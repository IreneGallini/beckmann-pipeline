import sys

from beckmann_nbo.parse_cmo import write_cn_ledger

if __name__ == "__main__":
    for mol_id in sys.argv[1:]:
        write_cn_ledger(mol_id)
