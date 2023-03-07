import os
import sys
import pickle
import pandas as pd
import collections
import gzip

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign


def confgenx_gen(work_path, optimizer):
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
            writer = Chem.SDWriter(f'{path}/input.sdf')
            mol = Chem.MolFromSmiles(cs)
            mol.SetProp('_Name', smiles)
            molH = Chem.AddHs(mol)
            writer.write(molH)
            writer.flush()
            writer.close()

            conf_number = number*2
            cmd = f'module load schrodinger && cd {path} && confgenx -n {conf_number} -m {conf_number} -WAIT input.sdf'
            os.system(cmd)
            suppl = Chem.MaeMolSupplier(gzip.open(f'{path}/input-out.maegz'), removeHs=False)
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
    confgenx_gen(work_path, optimizer)
