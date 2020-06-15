import subprocess
import argparse
import time
import os

import pandas as pd 
from rdkit import Chem
from rdkit.Chem import AllChem

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input xlsx file ", \
            required=True, type=str)

    args = parser.parse_args()
            
    filename = args.file
                    
    data = pd.read_excel(filename)

    molatomtypes = {}
    atomtypesset = set()

    dim = len(data["SMILES"])
    for idx, ss in enumerate(data["SMILES"]):
        start = time.time()

        mol = Chem.MolFromSmiles(str(ss))

        AllChem.Compute2DCoords(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d) 
        Chem.Kekulize(mol)
        mol_3D = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3D, randomSeed=0xf00d)

        basename = "molecule_"+str(idx+1)

        fout = Chem.SDWriter(basename+".mol")
        fout.write(mol_3D)
        fout.close()

        captured = ""
        logd = data[data["SMILES"] == ss]["LogD"].values[0]
        result = subprocess.run("obabel -imol "+  basename+ ".mol " +
            "-omol2 -O " + basename+".mol2", shell=True, check=True, \
            stdout=subprocess.PIPE, universal_newlines=True)
        #print(result.stdout)

        result = subprocess.run("./featuresbuild " +  basename+".mol2",  \
                shell=True, check=True, stdout=subprocess.PIPE, \
                universal_newlines=True) 
        atomtypes = result.stdout.split("\n")
        molatomtypes[str(ss)] = atomtypes
        for at in atomtypes:
            atomtypesset.add(at)

        os.remove(basename+".mol")
        os.remove(basename+".mol2")

        end = time.time()

        print("Mol %10d of %10d has %5d "%(idx+1, dim, mol.GetNumAtoms()), 
                " atoms and LogD %10.5f (%12.7s s)"%(logd, (end - start)))

