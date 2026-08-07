# Demo checklist: both pipelines, in person

## AIMNet2-PySCF

```bash
beckmann-pyscf predict --smiles "O=C1CCC2=CC=CC=C21" --name demo --plot
```

conformers → AIMNet2 optimization → PySCF scan, ~7 points at ~1-2 min each → prediction

When done, open the two result files it just wrote:
```bash
open beckmann_pyscf_runs/demo/wcnmax_vs_rno.png
cat beckmann_pyscf_runs/demo/summary.txt
```

## 3. Stage-by-stage
check if the scripts missed something:

```bash
beckmann-pyscf conformers --smiles "O=C1CCC2=CC=CC=C21" --name demo2
beckmann-pyscf optimize --conformers-sdf <path printed above>
```

Point out `optimize`'s printed output specifically: the AIMNet2 energy and
the resolved oxime atom map (C/N/O/aryl/alkyl indices). `scan`'s
`--ci`/`--ni`/`--oi`/`--c-aryl`/`--c-alkyl` flags 

## 4. Gaussian/NBO7 CLI 


```bash
beckmann-nbo init      # if not already done on her machine
beckmann-nbo verify    # confirms SSH reachable, g16 executable, NBO7 wrapper set up
beckmann-nbo predict --smiles "O=C1CCC2=C1C=CC=C2" --name demo3
```

Point out the printed `--dir` path and explain: copy that path into every
later command (`status`, `predict --continue`, `report`)

Then, without waiting on Citadel, show a finished benchmark molecule's
output:

```bash
beckmann-nbo --mol 002 status
beckmann-nbo --mol 002 report --out /tmp/demo_report_002 --advanced
open /tmp/demo_report_002/mol_002_E/wcnmax_vs_rno.png
```

(No `--dir` needed here: `002` is a benchmark molecule, so the CLI's
default directory already finds it.)

## 5. Where things live (~1 min)

Briefly note `INSTRUCTIONS.md` Section 4 (standalone scripts): these are
the literal scripts that built the 34-molecule benchmark set, most are
hardwired to that specific set (not general-purpose), and are not the
entry point for her own new molecules; the CLIs in Sections 1-2 are.
