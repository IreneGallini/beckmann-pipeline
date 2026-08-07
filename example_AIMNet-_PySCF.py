'''
Example running AIMNet-PySCF pipeline on a single SMILES string
'''
from beckmann_pyscf.pipeline import predict

if __name__ == "__main__":
    result = predict("O=C1CCC2=CC=CC=C21")  # example SMILES
    print(result["prediction"])