import pandas as pd 
import argparse
from rdkit import Chem

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input xlsx file ", \
            required=True, type=str)

    args = parser.parse_args()
            
    filename = args.file
                    
    data = pd.read_excel(filename)

    for ss in data["SMILES"]:
        mol = Chem.MolFromSmiles(str(ss))
        logd = data[data["SMILES"] == ss]["LogD"].values[0]
        print("Mol has %5d "%(mol.GetNumAtoms()), " atoms and LogD %10.5f"%(logd))

