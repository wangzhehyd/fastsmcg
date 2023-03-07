import os
import pickle
import random
import sys

import numpy as np
import pandas as pd
import torch
import yaml
from easydict import EasyDict
from model.featurization import featurize_mol_from_smiles
from model.inference import construct_conformers
from model.model import GeoMol
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from torch_geometric.data import Batch


def geomol_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    config_file = 'script/geomol/config.yml'

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    seed = 2022
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    model = GeoMol(**config)

    state_dict = torch.load(f'script/geomol/best_model.pt', map_location='cpu')
    model.load_state_dict(state_dict, strict=True)
    model.eval()

    input_lines = pd.read_csv(input_file, sep=',', skiprows=14, header=0).values
    mol_chunks = [input_lines[i:i + 100] for i in range(0, len(input_lines), 100)]

    writer = Chem.SDWriter(output_file_sdf)

    for chunk_index, chunk in enumerate(mol_chunks):
        pickle_slice = {}  # 类似这样 {'mol1': rdkit.Chem.rdchem.Mol, 'mol2': rdkit.Chem.rdchem.Mol}
        output_file_pkl = os.path.join(work_path, f'fastsmcg_output_slice{chunk_index:0>6}.pkl')
        for index, (smiles, title, num_conf) in enumerate(chunk):
            mol_title = f'{title}_{index + 1 + chunk_index * 100}'
            try:
                tg_data = featurize_mol_from_smiles(smiles, dataset='drugs')
                data = Batch.from_data_list([tg_data])
                model(data, inference=True, n_model_confs=num_conf)
                model_coords = construct_conformers(data, model)

                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                mol.SetProp('_Name', mol_title)
                n_atoms = tg_data.x.size(0)

                for x in model_coords.split(1, dim=1):
                    coords = x.squeeze(1).double().cpu().detach().numpy()
                    conf = Chem.Conformer(mol.GetNumAtoms())
                    for i in range(n_atoms):
                        conf.SetAtomPosition(i, Geometry.Point3D(coords[i, 0], coords[i, 1], coords[i, 2]))
                    mol.AddConformer(conf, assignId=True)
                rdMolAlign.AlignMolConformers(mol)

                for conf_id in range(mol.GetNumConformers()):
                    if optimizer == '1':
                        AllChem.UFFOptimizeMolecule(mol, maxIters=500, confId=conf_id)
                    elif optimizer == '2':
                        AllChem.MMFFOptimizeMolecule(mol, maxIters=500, confId=conf_id)
                    elif optimizer == '3':
                        AllChem.MMFFOptimizeMolecule(mol, maxIters=500, confId=conf_id, mmffVariant='MMFF94s')
                    writer.write(mol, confId=conf_id)
                    pickle_slice[mol_title] = mol
            except:
                pickle_slice[mol_title] = []

        with open(output_file_pkl, 'wb') as f:
            pickle.dump(pickle_slice, f)

    writer.flush()
    writer.close()


if __name__ == '__main__':
    work_path = sys.argv[1]
    optimizer = sys.argv[2]
    geomol_gen(work_path, optimizer)
