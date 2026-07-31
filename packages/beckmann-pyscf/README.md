# beckmann-pyscf

Open-source, HPC-free Beckmann rearrangement predictor: SMILES → AIMNet2 geometry → PySCF-native wCNmax scan → R/F prediction, with a Flask web UI. No Gaussian, no NBO7, no Citadel — verified by `backend/tests/test_no_hpc_dependency.py`, not just intended.

Depends on `beckmann-core` (pinned version) for the shared oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum R/F rule. `backend/beckmann_pyscf/engine/` is the validated PySCF wCNmax computation engine, relocated unchanged from the original prototype.

```bash
pip install -e .   # pulls in beckmann-core automatically
cd backend && python app.py   # http://localhost:5001

python -m pytest tests/ backend/tests/ -v -m "not slow" -c backend/pytest.ini
```

See the monorepo root `CLAUDE.md` for the full architecture writeup, including the validation caveat on how this package's wCNmax numbers compare to `research/pyscf_validation/`'s Gaussian-geometry comparison.
