# beckmann-core

Shared, method-agnostic library for the Beckmann rearrangement products (`beckmann-nbo`, `beckmann-pyscf`): ketone→oxime conversion, Auto3D conformer generation, AIMNet2/ASE geometry optimization, shared geometry primitives, the classical anti-periplanar baseline, and the wCNmax-minimum R/F rule.

No filesystem path conventions of its own — every function takes explicit `Path` arguments. See the monorepo root `CLAUDE.md` for the full architecture writeup.

```bash
pip install -e .
python -m pytest tests/ -v -m "not slow"
```
