# Descriptor summary (mol_001, mol_002, mol_003, mol_004, mol_005, mol_006, mol_007, mol_008, mol_009, mol_010, mol_011, mol_012, mol_013, mol_014, mol_015, mol_016, mol_017, mol_018, mol_019, mol_020, mol_021, mol_022, mol_023, mol_024, mol_025, mol_026, mol_027, mol_028, mol_029, mol_030, mol_031, mol_032, mol_033, mol_034)

d/dR = least-squares slope over each molecule's N-O scan series. 'wCNmax extremum' = interior local min/max in the wCNmax(R) series (the paper's central signature -- Table 2 reports this only for the one rearranging reference compound, none of the three fragmenting ones). 'crossing classification' = beckmann.dft.parse_cmo.classify_crossing()'s verdict on whether the wCNmax MO handoff is a CONFIRMED avoided crossing (bracketed narrow eigenvalue gap between the identity-tracked pre/post-handoff MO pair, AND roughly conserved CN weight across them) vs. an unconfirmed handoff vs. no handoff at all. Across all 34 molecules the crossing partner is consistently the N-O sigma*/sigma antibond, not the aryl-migrating C-C antibond (data/output/analysis/cn_crossing_report.csv), so no aryl-coefficient swap is checked.

| mol | exp | d(Ψ)/dR | d(log₁₀Λ)/dR | d(wCNmax)/dR | d(w17max)/dR | d(w78max)/dR | wCNmax extremum | dip depth | crossing classification |
|---|---|---|---|---|---|---|---|---|---|
| mol_001_E | F | 2.208 | 0.819 | 0.025 | -0.262 | 0.006 | yes @ R=1.6596 (MO 40, epsilon=-0.0099 a.u.) | 0.0052 | confirmed avoided crossing |
| mol_002_E | F | 2.018 | 0.429 | 0.021 | -0.174 | -0.026 | yes @ R=1.6592 (MO 49, epsilon=-0.0022 a.u.) | 0.0164 | confirmed avoided crossing |
| mol_003_E | F | 2.094 | -0.451 | 0.022 | 0.116 | -0.001 | yes @ R=1.6678 (MO 45, epsilon=-0.0094 a.u.) | 0.0078 | confirmed avoided crossing |
| mol_004_E | F | 2.153 | 0.646 | 0.019 | -0.286 | -0.033 | no | n/a | confirmed avoided crossing |
| mol_005_E | F | 2.154 | 0.485 | 0.019 | -0.160 | 0.019 | no | n/a | confirmed avoided crossing |
| mol_006_E | R | 2.066 | 0.287 | 0.021 | -0.104 | -0.021 | yes @ R=1.6608 (MO 45, epsilon=-0.0063 a.u.) | 0.1010 | confirmed avoided crossing |
| mol_007_E | F | 2.188 | -0.586 | 0.021 | 0.023 | -0.096 | no | n/a | confirmed avoided crossing |
| mol_008_E | F | 2.103 | -0.488 | 0.026 | 0.096 | -0.030 | yes @ R=1.5977 (MO 52, epsilon=0.0064 a.u.) | 0.0062 | confirmed avoided crossing |
| mol_009_E | R | 2.288 | -0.294 | 0.033 | 0.270 | 0.114 | yes @ R=1.6603 (MO 44, epsilon=-0.0089 a.u.) | 0.0192 | confirmed avoided crossing |
| mol_010_E | F | 2.219 | 0.493 | 0.022 | -0.272 | -0.091 | no | n/a | confirmed avoided crossing |
| mol_011_E | F | 2.154 | -0.937 | 0.022 | 0.190 | -0.129 | no | n/a | confirmed avoided crossing |
| mol_012_E | F | 2.089 | -0.529 | 0.014 | 0.091 | -0.049 | yes @ R=1.6989 (MO 53, epsilon=-0.0200 a.u.) | 0.0031 | confirmed avoided crossing |
| mol_013_E | R | 2.219 | 0.427 | 0.024 | 0.116 | 0.149 | yes @ R=1.6664 (MO 45, epsilon=-0.0060 a.u.) | 0.0065 | confirmed avoided crossing |
| mol_014_Z | F | 0.398 | -0.308 | -0.064 | 0.070 | 0.003 | yes @ R=1.7034 (MO 45, epsilon=-0.0138 a.u.) | 0.0439 | confirmed avoided crossing |
| mol_015_E | R | 2.218 | -0.784 | 0.021 | 0.317 | 0.045 | no | n/a | confirmed avoided crossing |
| mol_016_E | F | 2.211 | 0.829 | 0.098 | -0.187 | 0.055 | yes @ R=1.6573 (MO 48, epsilon=-0.0081 a.u.) | 0.1435 | confirmed avoided crossing |
| mol_017_E | F | 2.190 | 2.066 | 0.016 | -0.393 | 0.085 | no | n/a | confirmed avoided crossing |
| mol_018_E | F | 1.977 | -0.200 | 0.004 | 0.171 | 0.073 | yes @ R=1.6480 (MO 52, epsilon=-0.0111 a.u.) | 0.0331 | confirmed avoided crossing |
| mol_019_E | R | 1.973 | -0.392 | 0.047 | -0.076 | -0.094 | yes @ R=1.6505 (MO 44, epsilon=-0.0117 a.u.) | 0.0928 | confirmed avoided crossing |
| mol_020_E | R | 2.037 | -0.127 | 0.060 | 0.044 | -0.003 | yes @ R=1.6629 (MO 49, epsilon=-0.0016 a.u.) | 0.1212 | confirmed avoided crossing |
| mol_021_E | R | 1.749 | 0.072 | -0.012 | -0.110 | -0.085 | yes @ R=1.6533 (MO 52, epsilon=-0.0061 a.u.) | 0.2131 | confirmed avoided crossing |
| mol_022_E | R | 1.891 | -1.226 | 0.057 | 0.069 | -0.230 | yes @ R=1.6497 (MO 48, epsilon=-0.0112 a.u.) | 0.1116 | confirmed avoided crossing |
| mol_023_E | R | 2.005 | -1.241 | 0.110 | 0.269 | 0.017 | yes @ R=1.5940 (MO 56, epsilon=0.0071 a.u.) | 0.0198 | confirmed avoided crossing |
| mol_024_E | R | 2.022 | 0.006 | 0.124 | -0.041 | -0.031 | yes @ R=1.6537 (MO 48, epsilon=-0.0119 a.u.) | 0.1393 | confirmed avoided crossing |
| mol_025_E | R | 2.060 | -0.015 | 0.048 | -0.100 | -0.100 | yes @ R=1.6512 (MO 52, epsilon=-0.0129 a.u.) | 0.0950 | confirmed avoided crossing |
| mol_026_E | R | 2.003 | 0.479 | 0.008 | -0.171 | -0.076 | yes @ R=1.6492 (MO 48, epsilon=-0.0164 a.u.) | 0.0529 | confirmed avoided crossing |
| mol_027_E | R | 1.969 | -1.223 | -0.034 | 0.122 | -0.082 | yes @ R=1.6994 (MO 57, epsilon=-0.0131 a.u.) | 0.0517 | confirmed avoided crossing |
| mol_028_E | R | 2.225 | -0.865 | 0.396 | 0.259 | -0.028 | yes @ R=1.6668 (MO 49, epsilon=0.0044 a.u.) | 0.0533 | confirmed avoided crossing |
| mol_029_Z | R | 0.407 | -0.747 | 0.024 | 0.079 | -0.065 | yes @ R=1.6566 (MO 49, epsilon=-0.0036 a.u.) | 0.2118 | confirmed avoided crossing |
| mol_030_E | R | 2.011 | 0.095 | 0.206 | 0.015 | 0.025 | yes @ R=1.6564 (MO 52, epsilon=-0.0132 a.u.) | 0.1765 | confirmed avoided crossing |
| mol_031_E | R | 1.915 | 0.087 | 0.125 | -0.200 | -0.107 | yes @ R=1.6511 (MO 48, epsilon=-0.0153 a.u.) | 0.0962 | confirmed avoided crossing |
| mol_032_E | R | 1.924 | -1.309 | -0.036 | 0.107 | -0.254 | yes @ R=1.6982 (MO 62, epsilon=-0.0117 a.u.) | 0.0330 | confirmed avoided crossing |
| mol_033_E | R | 1.902 | 0.171 | -0.019 | 0.006 | 0.048 | yes @ R=1.6487 (MO 52, epsilon=-0.0157 a.u.) | 0.0377 | confirmed avoided crossing |
| mol_034_E | R | -0.051 | -0.863 | 0.186 | 0.046 | -0.128 | yes @ R=1.6690 (MO 49, epsilon=-0.0001 a.u.) | 0.0678 | confirmed avoided crossing |
