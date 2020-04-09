import pandas as pd 
import argparse
import pybel

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="input xlsx file ", \
            required=True, type=str)

    args = parser.parse_args()
            
    filename = args.file
                    
    data = pd.read_excel(filename, sheet_name="S1")

    for ss in data["SMILES structure"]:
        mol = pybel.readstring("smi", ss)

