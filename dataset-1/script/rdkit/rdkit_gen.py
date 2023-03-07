import os
import sys
import pickle
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign


def rdkit_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    input_lines = pd.read_csv(input_file, sep=',', skiprows=14, header=0).values

    writer = Chem.SDWriter(output_file_sdf)
    params = rdDistGeom.ETKDGv3()
    params.useSmallRingTorsions = True
    params.pruneRmsThresh = 0.0

    mol_checks = [input_lines[i:i + 100] for i in range(0, len(input_lines), 100)]
    for check_index, check in enumerate(mol_checks):
        pickle_slice = {}  # 类似这样 {'mol1': rdkit.Chem.rdchem.Mol, 'mol2': rdkit.Chem.rdchem.Mol}
        output_file_pkl = os.path.join(work_path, f'fastsmcg_output_slice{check_index:0>6}.pkl')

        for index, (smiles, title, number) in enumerate(check):
            mol_title = f'{title}_{index + 1 + check_index * 100}'
            try:
                mol = Chem.MolFromSmiles(smiles)
                mol.SetProp('_Name', mol_title)
                molH = Chem.AddHs(mol)
                conf_ids = rdDistGeom.EmbedMultipleConfs(molH, number, params)
                rdMolAlign.AlignMolConformers(molH)
                for conf_id in conf_ids:
                    if optimizer == '1':
                        AllChem.UFFOptimizeMolecule(molH, maxIters=500, confId=conf_id)
                    elif optimizer == '2':
                        AllChem.MMFFOptimizeMolecule(molH, maxIters=500, confId=conf_id)
                    elif optimizer == '3':
                        AllChem.MMFFOptimizeMolecule(mol, maxIters=500, confId=conf_id, mmffVariant='MMFF94s')
                    writer.write(molH, confId=conf_id)
                pickle_slice[mol_title] = molH
            except:
                pickle_slice[mol_title] = []

        with open(output_file_pkl, 'wb') as f:
            pickle.dump(pickle_slice, f)

    writer.flush()
    writer.close()

if __name__ == '__main__':
    work_path = sys.argv[1]
    optimizer = sys.argv[2]
    rdkit_gen(work_path, optimizer)
