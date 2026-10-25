"""
Builds the frontend as a plain static site for GitHub Pages
(.github/workflows/deploy-pages.yml), so that the page and its 24 MB
Ketcher bundle are served for free by GitHub instead of by Cloud Run,
whose outbound data is billed past 1 GiB a month. The page then calls
the Cloud Run API cross-origin (Cloud Run allows it via
BECKMANN_CORS_ORIGINS).

    python build_static.py --api-base https://beckmann-pyscf-xxxx-uc.a.run.app --out _site

Renders templates/index.html (the same template Flask serves locally) in
sync mode with the given API base, and copies static/ next to it.
Needs only jinja2.
"""
import argparse
import shutil
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

HERE = Path(__file__).resolve().parent


def build(api_base: str, out: Path) -> None:
    api_base = api_base.rstrip("/")
    local = api_base.startswith(("http://localhost", "http://127.0.0.1"))  # local testing
    if not (api_base.startswith("https://") or local):
        raise SystemExit(f"--api-base must be an https:// URL, got {api_base!r}")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    env = Environment(loader=FileSystemLoader(HERE / "templates"), autoescape=True)
    html = env.get_template("index.html").render(api_mode="sync", api_base=api_base)
    (out / "index.html").write_text(html)
    shutil.copytree(HERE / "static", out / "static")
    # Serve files as-is: no Jekyll processing on GitHub Pages.
    (out / ".nojekyll").write_text("")
    print(f"built {out} (API: {api_base})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--api-base", required=True, help="Cloud Run service URL")
    parser.add_argument("--out", type=Path, default=HERE / "_site")
    args = parser.parse_args()
    build(args.api_base, args.out)
