"""
Job store + executor for the async submit/poll pattern. Each job runs in
its own OS process (multiprocessing, spawn context), not a thread -- only a
process can be forcibly reclaimed on timeout, which is the whole point:
AIMNet2/PySCF computations run on real CPU for real minutes, and a stuck
one has to actually stop, not just get reported as stopped while it keeps
burning CPU in the background.

Two roles, split so any number of gunicorn workers can serve HTTP:

- **Web side** (every gunicorn worker): `submit()` and `get_status()`. These
  only ever touch the job *store* -- never a Process object.
- **Dispatcher** (exactly one process): `_reap_once()`, looped. Owns every
  job's multiprocessing.Process (_PROCESSES/_STARTED_AT below are local to
  it), promotes pending jobs into free slots, collects results, enforces
  timeouts.

The store is picked by `settings.REDIS_URL`:

- unset -> `MemoryStore`: plain in-process dicts. The dispatcher is the
  reaper thread `start_reaper()` runs inside the one web process, so this
  mode still requires `gunicorn -w 1` (local dev, `python app.py`, tests).
- set   -> `RedisStore`: job records/queue live in Redis, the dispatcher is
  its own process (`python dispatcher.py`), and gunicorn can run as many
  workers as it likes.

Statuses: "queued" | "running" | "done" | "error" | "invalid_input" | "timeout".
"invalid_input" and "timeout" are distinct from a generic "error" on
purpose -- a bad SMILES and a stuck SCF convergence are different problems
and should read differently to the user.
"""
import importlib
import json
import logging
import multiprocessing
import os
import queue as queue_module
import signal
import sys
import threading
import time
import uuid

from beckmann_pyscf import settings

logger = logging.getLogger("beckmann_pyscf.jobs")

_ctx = multiprocessing.get_context("spawn")

_TERMINAL_STATUSES = {"done", "error", "invalid_input", "timeout"}


# --------------------------------------------------------------------------
# Stores. Both expose the same small interface; neither ever holds a
# Process/Queue or anything else non-JSON-serializable -- a job record is
# exactly what /status exposes via jsonify.
# --------------------------------------------------------------------------

class MemoryStore:
    """Single-process store: plain dicts guarded by one lock."""

    def __init__(self):
        self._lock = threading.Lock()
        self._jobs: dict[str, dict] = {}
        self._pending: list[dict] = []
        self._running: set[str] = set()
        self._finished_at: dict[str, float] = {}

    def try_enqueue(self, job_id: str, entry: dict, capacity: int) -> bool:
        with self._lock:
            if len(self._pending) + len(self._running) >= capacity:
                return False
            self._jobs[job_id] = {"status": "queued"}
            self._pending.append(entry)
            return True

    def pop_pending(self) -> dict | None:
        with self._lock:
            return self._pending.pop(0) if self._pending else None

    def mark_running(self, job_id: str) -> None:
        with self._lock:
            self._running.add(job_id)
            self._jobs[job_id] = {"status": "running"}

    def finish(self, job_id: str, record: dict) -> None:
        with self._lock:
            self._running.discard(job_id)
            self._jobs[job_id] = record
            self._finished_at[job_id] = time.time()

    def get(self, job_id: str) -> dict | None:
        with self._lock:
            return self._jobs.get(job_id)

    def running_ids(self) -> list[str]:
        with self._lock:
            return list(self._running)

    def job_count(self) -> int:
        with self._lock:
            return len(self._jobs)

    def sweep(self) -> None:
        now = time.time()
        with self._lock:
            expired = [
                job_id for job_id, finished in self._finished_at.items()
                if now - finished > settings.JOB_RETENTION_SECONDS
            ]
            for job_id in expired:
                self._jobs.pop(job_id, None)
                self._finished_at.pop(job_id, None)

    def ping(self) -> bool:
        return True


class RedisStore:
    """Shared store: job records are `job:<id>` JSON strings (terminal ones
    carry a TTL of JOB_RETENTION_SECONDS, so Redis does the sweeping),
    pending jobs a FIFO list, running job ids a set."""

    def __init__(self, client, prefix: str = "beckmann:"):
        self.r = client
        self._pending_key = prefix + "pending"
        self._running_key = prefix + "running"
        self._job_prefix = prefix + "job:"

    def _job_key(self, job_id: str) -> str:
        return self._job_prefix + job_id

    def try_enqueue(self, job_id: str, entry: dict, capacity: int) -> bool:
        # WATCH/MULTI so the capacity check and the push are atomic: two
        # gunicorn workers racing /submit can't both squeeze into the last
        # slot. A concurrent change to either watched key aborts the
        # transaction and transaction() retries it.
        def _txn(pipe):
            if pipe.llen(self._pending_key) + pipe.scard(self._running_key) >= capacity:
                pipe.multi()
                return False
            pipe.multi()
            pipe.set(self._job_key(job_id), json.dumps({"status": "queued"}))
            pipe.rpush(self._pending_key, json.dumps(entry))
            return True

        return self.r.transaction(
            _txn, self._pending_key, self._running_key, value_from_callable=True
        )

    def pop_pending(self) -> dict | None:
        raw = self.r.lpop(self._pending_key)
        return json.loads(raw) if raw is not None else None

    def mark_running(self, job_id: str) -> None:
        pipe = self.r.pipeline()
        pipe.sadd(self._running_key, job_id)
        pipe.set(self._job_key(job_id), json.dumps({"status": "running"}))
        pipe.execute()

    def finish(self, job_id: str, record: dict) -> None:
        pipe = self.r.pipeline()
        pipe.srem(self._running_key, job_id)
        pipe.set(self._job_key(job_id), json.dumps(record),
                 ex=max(1, settings.JOB_RETENTION_SECONDS))
        pipe.execute()

    def get(self, job_id: str) -> dict | None:
        raw = self.r.get(self._job_key(job_id))
        return json.loads(raw) if raw is not None else None

    def running_ids(self) -> list[str]:
        return list(self.r.smembers(self._running_key))

    def job_count(self) -> int:
        return sum(1 for _ in self.r.scan_iter(match=self._job_prefix + "*"))

    def sweep(self) -> None:
        pass  # terminal records expire via their own TTL

    def ping(self) -> bool:
        try:
            return bool(self.r.ping())
        except Exception:
            return False


_store = None


def get_store():
    """The process-wide store, created on first use from settings.REDIS_URL."""
    global _store
    if _store is None:
        if settings.REDIS_URL:
            import redis
            _store = RedisStore(redis.Redis.from_url(settings.REDIS_URL, decode_responses=True))
        else:
            _store = MemoryStore()
    return _store


def uses_shared_store() -> bool:
    return isinstance(get_store(), RedisStore)


# --------------------------------------------------------------------------
# Dispatcher-local process bookkeeping. Only the dispatcher process ever
# populates these.
# --------------------------------------------------------------------------

_PROCESSES: dict[str, multiprocessing.Process] = {}
_STARTED_AT: dict[str, float] = {}
_RESULT_QUEUE = _ctx.Queue()
_LOCK = threading.Lock()

_reaper_started = False


def _worker(job_id: str, smiles: str, result_queue) -> None:
    """Runs IN the child process -- imports beckmann_pyscf.pipeline here
    (not at this module's top level) so merely importing jobs.py in the
    parent doesn't eagerly pull in torch/rdkit/pyscf before a job is even
    submitted. Never touches the store directly; the only way results get
    back to the dispatcher is this queue."""
    from beckmann_pyscf.pipeline import predict

    try:
        result = predict(smiles)
        result_queue.put({"job_id": job_id, "status": "done", "result": result})
    except ValueError as e:
        result_queue.put({"job_id": job_id, "status": "invalid_input", "error": str(e)})
    except Exception as e:
        result_queue.put({"job_id": job_id, "status": "error", "error": str(e)})


def _run_in_own_group(target, job_id: str, smiles: str, result_queue) -> None:
    """The actual multiprocessing entrypoint -- wraps `target` so every job
    process becomes the leader of its own new process group (os.setpgrp())
    before anything else runs. A real smoke test found that AIMNet2/torch
    spawn their own internal multiprocessing workers as grandchildren of our
    tracked process -- os.kill()/proc.terminate() aimed at just that one pid
    doesn't reach them, so they survive a "successful" timeout kill. Giving
    the job its own pgid lets the reaper target the whole group at once
    (os.killpg) instead."""
    os.setpgrp()
    target(job_id, smiles, result_queue)


def _target_ref(target) -> str:
    """A pending job's target has to cross the store as a string (it may
    sit in Redis until a different process -- the dispatcher -- starts it)."""
    return f"{target.__module__}:{target.__qualname__}"


def _resolve_target(ref: str):
    module_name, _, name = ref.partition(":")
    module = sys.modules.get(module_name) or importlib.import_module(module_name)
    return getattr(module, name)


def _start(job_id: str, smiles: str, target) -> None:
    """Caller must hold _LOCK. `target` is normally the real _worker (which
    runs the actual predict() pipeline); tests override it with a fast
    dummy so the job-store mechanics (queueing, timeout, expiry) can be
    exercised in milliseconds instead of the real multi-minute pipeline.
    Must be a top-level, picklable function either way -- multiprocessing's
    spawn context re-imports it by (module, name) in the child. daemon=False:
    the real predict() pipeline (Auto3D/AIMNet2) spawns its own child
    processes internally, and Python forbids a daemonic process from having
    children -- confirmed by a real "daemonic processes are not allowed to
    have children" job failure during manual smoke testing. Non-daemon means
    an abnormal (SIGKILL) death of the parent won't auto-reap these, but
    normal shutdown (SIGTERM) still lets them be terminated explicitly,
    which _kill_group(), timeout handling and run_dispatcher() already do."""
    process = _ctx.Process(
        target=_run_in_own_group, args=(target, job_id, smiles, _RESULT_QUEUE), daemon=False
    )
    process.start()
    _PROCESSES[job_id] = process
    _STARTED_AT[job_id] = time.time()
    get_store().mark_running(job_id)
    logger.info("job %s started (pid=%s)", job_id, process.pid)


def _promote() -> None:
    """Start pending jobs while there are free slots. Dispatcher only."""
    store = get_store()
    with _LOCK:
        while len(_PROCESSES) < settings.MAX_CONCURRENT_JOBS:
            entry = store.pop_pending()
            if entry is None:
                break
            try:
                target = _resolve_target(entry["target"])
            except Exception as e:
                logger.exception("job %s: cannot resolve target %r", entry["job_id"], entry["target"])
                store.finish(entry["job_id"], {"status": "error", "error": f"internal error: {e}"})
                continue
            _start(entry["job_id"], entry["smiles"], target)


def submit(smiles: str, target=_worker) -> str:
    """Queues a background prediction job, returns its job_id. Raises
    RuntimeError if the server is already at capacity (MAX_CONCURRENT_JOBS
    running + MAX_QUEUE_DEPTH waiting) -- app.py turns that into a 503.
    `target` is a testing hook -- see _start().

    With the in-process MemoryStore this process *is* the dispatcher, so a
    free slot is filled immediately rather than on the reaper's next pass."""
    job_id = str(uuid.uuid4())
    store = get_store()
    entry = {"job_id": job_id, "smiles": smiles, "target": _target_ref(target)}
    capacity = settings.MAX_CONCURRENT_JOBS + settings.MAX_QUEUE_DEPTH
    if not store.try_enqueue(job_id, entry, capacity):
        raise RuntimeError("Server is at capacity, try again shortly")
    if not uses_shared_store():
        _promote()
    if store.get(job_id) == {"status": "queued"}:
        logger.info("job %s queued", job_id)
    return job_id


def get_status(job_id: str) -> dict | None:
    """The job's current {"status", ...} dict, or None if job_id is unknown
    (never submitted, or already expired past JOB_RETENTION_SECONDS)."""
    return get_store().get(job_id)


def _kill_group(job_id: str, proc) -> None:
    # Target the whole process group (pid == pgid: _run_in_own_group
    # made this job its own group leader via os.setpgrp()), not just
    # the single tracked pid. A real smoke test found AIMNet2/torch's
    # own internal multiprocessing workers -- grandchildren of our
    # tracked process -- survive a plain os.kill()/proc.terminate()
    # aimed at just the direct child. Also don't trust
    # proc.is_alive()'s word before escalating: the same smoke test
    # saw it report False (skipping straight past a kill attempt)
    # while the OS process group was still alive minutes later, a
    # multiprocessing/spawn-context bookkeeping quirk on macOS.
    #
    # If the child hasn't reached its os.setpgrp() yet (killed within its
    # first moments -- a spawn child spends a while re-importing modules
    # first), there is no process group `pid` and killpg raises
    # ProcessLookupError while the process itself is alive and well; fall
    # back to signalling the pid directly then. It can't have spawned
    # grandchildren before setpgrp, so the pid alone is the whole job.
    pid = proc.pid

    def _signal(sig):
        try:
            os.killpg(pid, sig)
        except ProcessLookupError:
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
        except OSError:
            logger.exception("job %s: os.killpg(%s, %s) failed", job_id, pid, sig)

    if pid is not None:
        _signal(signal.SIGTERM)
    proc.join(timeout=5)
    if pid is not None:
        _signal(signal.SIGKILL)
    proc.join(timeout=5)
    if proc.is_alive():
        logger.error("job %s cleanup FAILED to terminate pid %s", job_id, pid)


def _reap_once() -> None:
    """One dispatcher pass: drain finished results, kill timed-out
    processes, promote queued jobs into freed slots, sweep expired
    terminal records."""
    store = get_store()

    while True:
        try:
            msg = _RESULT_QUEUE.get_nowait()
        except queue_module.Empty:
            break
        job_id = msg.pop("job_id")
        with _LOCK:
            proc = _PROCESSES.pop(job_id, None)
            _STARTED_AT.pop(job_id, None)
        store.finish(job_id, msg)
        if proc is not None:
            proc.join(timeout=5)
        logger.info("job %s finished: %s", job_id, msg["status"])

    now = time.time()
    with _LOCK:
        timed_out = [
            job_id for job_id, started in _STARTED_AT.items()
            if now - started > settings.JOB_TIMEOUT_SECONDS
        ]
    for job_id in timed_out:
        with _LOCK:
            proc = _PROCESSES.pop(job_id, None)
            _STARTED_AT.pop(job_id, None)
        if proc is not None:
            _kill_group(job_id, proc)
        store.finish(job_id, {"status": "timeout"})
        logger.warning("job %s timed out after %ss", job_id, settings.JOB_TIMEOUT_SECONDS)

    _promote()
    store.sweep()


def recover_orphans() -> None:
    """Dispatcher startup: any job the store still lists as running isn't
    actually running anywhere (its dispatcher died or the container was
    restarted) -- report it as failed instead of leaving clients polling
    "running" forever."""
    store = get_store()
    with _LOCK:
        orphans = [j for j in store.running_ids() if j not in _PROCESSES]
    for job_id in orphans:
        store.finish(job_id, {"status": "error",
                              "error": "Server restarted while this job was running, please resubmit"})
        logger.warning("job %s marked as error: orphaned by a dispatcher restart", job_id)


def _reaper_loop(poll_interval: float) -> None:
    while True:
        try:
            _reap_once()
        except Exception:
            logger.exception("reaper pass failed")
        time.sleep(poll_interval)


def start_reaper(poll_interval: float = 2.0) -> None:
    """Starts the in-process dispatcher thread -- MemoryStore mode only
    (with Redis, the dispatcher is its own process: run_dispatcher()).
    Idempotent -- app.py calls this once at startup; must NOT happen as a
    module-level import side effect, since multiprocessing (spawn)
    re-executes this module's top level in every child worker process too,
    which would otherwise spin up a redundant, useless reaper in each one."""
    global _reaper_started
    if _reaper_started or uses_shared_store():
        return
    _reaper_started = True
    threading.Thread(target=_reaper_loop, args=(poll_interval,), daemon=True).start()


def run_dispatcher(poll_interval: float = 1.0) -> None:
    """Foreground dispatcher loop for RedisStore mode (see dispatcher.py).
    On SIGTERM/SIGINT, kills every running job's process group before
    exiting, so a container stop doesn't leave orphaned PySCF processes."""
    if not uses_shared_store():
        raise SystemExit("run_dispatcher() requires BECKMANN_REDIS_URL to be set")

    def _shutdown(signum, frame):
        logger.info("dispatcher received signal %s, terminating %d running job(s)",
                    signum, len(_PROCESSES))
        with _LOCK:
            running = list(_PROCESSES.items())
            _PROCESSES.clear()
            _STARTED_AT.clear()
        for job_id, proc in running:
            _kill_group(job_id, proc)
            get_store().finish(job_id, {"status": "error",
                                        "error": "Server restarted while this job was running, please resubmit"})
        sys.exit(0)

    signal.signal(signal.SIGTERM, _shutdown)
    signal.signal(signal.SIGINT, _shutdown)

    recover_orphans()
    logger.info("dispatcher started (max_concurrent=%s, timeout=%ss)",
                settings.MAX_CONCURRENT_JOBS, settings.JOB_TIMEOUT_SECONDS)
    _reaper_loop(poll_interval)


def run_blocking(smiles: str, heartbeat_s: float = 20.0, target=_worker):
    """Synchronous mode (settings.API_MODE == "sync", used on Cloud Run):
    runs one job in its own process group *inside the caller's request*,
    bypassing the store/dispatcher entirely. A generator: yields None every
    `heartbeat_s` while the job runs, then yields the final record (same
    shape as a /status record: done/error/invalid_input/timeout) and stops.

    Closing the generator early -- which is what gunicorn/werkzeug do when
    the client disconnects and a heartbeat write fails -- kills the job's
    process group in the `finally` below, so an abandoned scan never keeps
    running (on Cloud Run it would sit CPU-throttled, holding the instance's
    one scan slot, until JOB_TIMEOUT_SECONDS)."""
    result_queue = _ctx.Queue()
    proc = _ctx.Process(
        target=_run_in_own_group, args=(target, "sync", smiles, result_queue), daemon=False
    )
    proc.start()
    logger.info("sync job started (pid=%s)", proc.pid)
    deadline = time.time() + settings.JOB_TIMEOUT_SECONDS
    finished = False
    try:
        while True:
            wait = min(heartbeat_s, max(0.0, deadline - time.time()))
            try:
                msg = result_queue.get(timeout=wait)
            except queue_module.Empty:
                if time.time() >= deadline:
                    _kill_group("sync", proc)
                    finished = True
                    logger.warning("sync job timed out after %ss", settings.JOB_TIMEOUT_SECONDS)
                    yield {"status": "timeout"}
                    return
                yield None
                continue
            msg.pop("job_id", None)
            proc.join(timeout=5)
            finished = True
            logger.info("sync job finished: %s", msg["status"])
            yield msg
            return
    finally:
        if not finished:
            logger.warning("sync job abandoned by client, killing pid %s", proc.pid)
            _kill_group("sync", proc)


def reset_for_tests(store=None) -> None:
    """Tests only: drop all state and install `store` (default: a fresh
    MemoryStore). Terminates anything a previous test left running."""
    global _store
    with _LOCK:
        leftovers = list(_PROCESSES.values())
        _PROCESSES.clear()
        _STARTED_AT.clear()
    for proc in leftovers:
        if proc.is_alive():
            proc.terminate()
            proc.join(timeout=5)
    while True:
        try:
            _RESULT_QUEUE.get_nowait()
        except queue_module.Empty:
            break
    _store = store if store is not None else MemoryStore()
