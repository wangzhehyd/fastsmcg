import os
import sys
import pickle
import pandas as pd
import collections

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign


def conformator_gen(work_path, optimizer):
    input_file = os.path.join(work_path, '../dataset-2.csv')
    output_file = os.path.join(work_path, 'output.pkl')

    input_lines = pd.read_csv(input_file, sep=',', header=0).values

    pickle_mols = collections.OrderedDict()
    for index, (smiles, number, cs) in enumerate(input_lines):
        mols = []
        path = f'{index:0>4}'
        try:
            cmd = f'mkdir {path}'
            os.system(cmd)
            with open(f'{path}/input.smi','w') as f:
                f.write(f'{cs}')
            cmd = f'conformator -i {path}/input.smi -o {path}/output.sdf -n {number*2} -q 1'
            os.system(cmd)
            suppl = Chem.SDMolSupplier(f'{path}/output.sdf', removeHs=False)
            for m in suppl:
                mols.append(m)
            pickle_mols[cs] = mols
        except:
            pickle_mols[cs] = []
    with open(output_file, 'wb') as f:
        pickle.dump(pickle_mols, f)

if __name__ == '__main__':
    work_path = sys.argv[1]
    optimizer = sys.argv[2]
    conformator_gen(work_path, optimizer)
