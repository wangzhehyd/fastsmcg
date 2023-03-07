import os
import pickle
import sys

import pandas as pd
import torch
import yaml
from easydict import EasyDict
from rdkit.Chem import rdMolAlign

from diffusion.sampling import *
from utils.utils import get_model


def embed_func(mol, numConfs):
    AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, numThreads=1)
    return mol


def sample_confs(smiles, n_confs, model, config):
    mol, data = get_seed(smiles, dataset='drugs')
    if not mol:
        return None

    n_rotable_bonds = int(data.edge_mask.sum())
    conformers, pdb = embed_seeds(mol, data, n_confs, single_conf=False, pdb=None, embed_func=embed_func, mmff=False)

    if not conformers:
        return None

    if n_rotable_bonds > 0.5:
        conformers = perturb_seeds(conformers, pdb)
        conformers = sample(conformers, model, config.sigma_max, config.sigma_min, config.inference_steps, config.batch_size, config.ode, config.likelihood, pdb)

    mols = [pyg_to_mol(mol, conf, None, rmsd=False) for conf in conformers]

    return mols


def torsional_diffusion_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    config_file = 'script/torsional_diffusion/config.yml'

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    batch_size = config.batch_size

    model = get_model(config)
    state_dict = torch.load('script/torsional_diffusion/best_model.pt', map_location=torch.device('cpu'))
    model.load_state_dict(state_dict, strict=True)
    model = model.to(device)
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
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                mol.SetProp('_Name', mol_title)
                for m in sample_confs(smiles, num_conf, model, config):
                    mol.AddConformer(m.GetConformer(0), assignId=True)

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
    torsional_diffusion_gen(work_path, optimizer)
