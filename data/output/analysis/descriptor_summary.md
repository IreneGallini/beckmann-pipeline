# Descriptor summary (mol_002, 006, 021)

mol_020 omitted from this run -- its relaxed scan hit the step cap on scan point 5 without converging (oscillating ring pucker) and the job crashed after the NBO7 single-point; see mol_020_E_scan.log. Will be added back once that scan is fixed and rerun.

d/dR = least-squares slope over the 5-point N-O scan. 'wCNmax extremum' = interior local min/max in the wCNmax(R) series (the paper's central signature -- Table 2 reports this only for the one rearranging reference compound, none of the three fragmenting ones).

| mol | exp | d(Ψ)/dR | d(log₁₀Λ)/dR | d(wCNmax)/dR | d(w17max)/dR | d(w78max)/dR | wCNmax extremum |
|---|---|---|---|---|---|---|---|
| mol_002_E | F | 2.337 | 0.673 | 0.020 | 0.027 | 0.154 | no |
| mol_006_E | R | 2.426 | -0.653 | 0.018 | n/a | -0.113 | no |
| mol_021_E | R | 2.001 | -0.832 | 0.027 | n/a | -0.147 | yes @ R=1.7033 (MO 53, epsilon=-0.0016 a.u.) |
