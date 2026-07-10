# Descriptor summary (mol_002, 006, 021)

mol_020 omitted from this run -- its relaxed scan hit the step cap on scan point 5 without converging (oscillating ring pucker) and the job crashed after the NBO7 single-point; see JOB_ISSUES.md in data/output/dft_opt/. Will be added back once that scan is fixed and rerun.

Lambda/log_lambda are left undefined (n/a) rather than floored when w17max isn't found anywhere in the LUMO..LUMO+0.4 a.u. virtual-orbital window (mol_006 and mol_021, all stages) -- an earlier version of this pipeline substituted a 1e-3 floor for the missing w17max, which inflated Lambda into a division-floor artifact (values of 50-100) rather than a real ratio. See beckmann/dft/parse_cmo.py docstring.

d/dR = least-squares slope over the 5-point N-O scan. 'wCNmax extremum' = interior local min/max in the wCNmax(R) series (the paper's central signature -- Table 2 reports this only for the one rearranging reference compound, none of the three fragmenting ones).

| mol | exp | d(Ψ)/dR | d(log₁₀Λ)/dR | d(wCNmax)/dR | d(w17max)/dR | d(w78max)/dR | wCNmax extremum |
|---|---|---|---|---|---|---|---|
| mol_002_E | F | 2.337 | 0.673 | 0.020 | 0.027 | 0.154 | no |
| mol_006_E | R | 2.426 | n/a | 0.018 | n/a | -0.113 | no |
| mol_021_E | R | 2.001 | n/a | 0.027 | n/a | -0.147 | yes @ R=1.7033 (MO 53, epsilon=-0.0016 a.u.) |
