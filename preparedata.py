import subprocess
import argparse
import time
import os

import numpy as np
import pandas as pd 
from rdkit import Chem
from rdkit.Chem import AllChem

import mllogdpcommon

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input xlsx file ", \
            required=True, type=str)
    parser.add_argument("-o","--outfilename", help="output filename default=descriptors.txt ", \
            required=False, type=str, default="descriptors.txt")

    args = parser.parse_args()
            
    filename = args.file
                    
    data = pd.read_excel(filename)

    molatomtypes = {}
    atomtypesset = set()
    mollogd = {}

    errorscounter = 0
    errorssmiles = []
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

        try:
            result = subprocess.run("./featuresbuild " +  basename+".mol2",  \
                shell=True, check=True, stdout=subprocess.PIPE, \
                universal_newlines=True) 

            atomtypes = result.stdout.split("\n")
            molatomtypes[str(ss)] = atomtypes
            mollogd[str(ss)] = logd
            for at in atomtypes:
                atomtypesset.add(at)
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, " for ", str(ss))
            errorssmiles.append(str(ss))
            errorscounter += 1

        os.remove(basename+".mol")
        os.remove(basename+".mol2")

        end = time.time()

        print("Mol %10d of %10d has %5d "%(idx+1, dim, mol.GetNumAtoms()), 
                " atoms and LogD %10.5f (%12.7s s)"%(logd, (end - start)))

    print(errorscounter, " errors out of ", dim, " molecules ")
    for i, s in enumerate(errorssmiles):
        print("  %5d "%(i+1), s)

    orderedset = []
    for e in atomtypesset:
        orderedset.append(e)
     
    fp = open(args.outfilename, "w")
    for mol in molatomtypes:
        atomtypes = molatomtypes[mol]
        desc = to_descriptor (orderedset, atomtypes)
        line = mol + " , "
        for v in desc:
            line += "%5d , "%(v)
        line += " %15.6f \n"%(mollogd[mol])
        fp.write(line)

    fp.close()

