# beckmann-nbo

Gaussian/NBO7/Citadel product for the Beckmann rearrangement pipeline: DFT input generation, HPC submission, NBO7 log parsing, and NBO-specific descriptors (Psi, Lambda). Depends on `beckmann-core` (pinned version) for the shared oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum R/F rule.

See the monorepo root `CLAUDE.md` for the full DFT/HPC workflow, NBO7 setup on Citadel, and output file reference.

```bash
pip install -e .   # pulls in beckmann-core automatically
python -m pytest tests/ -v
```

Config (`.env`, Citadel connection details) lives at the monorepo root — copy `.env.example` there and fill in values.
