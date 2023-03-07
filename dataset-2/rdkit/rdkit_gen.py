import os
import sys
import pickle
import pandas as pd
import collections

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign


def rdkit_gen(work_path, optimizer):
    input_file = os.path.join(work_path, '../dataset-2.csv')
    output_file = os.path.join(work_path, 'output.pkl')

    input_lines = pd.read_csv(input_file, sep=',', header=0).values

    params = rdDistGeom.ETKDGv3()
    params.useSmallRingTorsions = True
    params.numThreads = 0
    #params.pruneRmsThresh = 0.5

    pickle_mols = collections.OrderedDict()
    for index, (smiles, number, cs) in enumerate(input_lines):
        mols = []
        path = f'{index:0>4}'
        try:
            cmd = f'mkdir {path}'
            os.system(cmd)
            writer = Chem.SDWriter(f'{path}/output.sdf')

            mol = Chem.MolFromSmiles(cs)
            mol.SetProp('_Name', smiles)
            molH = Chem.AddHs(mol)
            conf_ids = rdDistGeom.EmbedMultipleConfs(molH, number*2, params)
            rdMolAlign.AlignMolConformers(molH)

            for conf_id in conf_ids:
                writer.write(molH, confId=conf_id)

            writer.flush()
            writer.close()
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
    rdkit_gen(work_path, optimizer)
