import copy
import os
import pickle
import random
import sys

import numpy as np
import pandas as pd
import torch
import torch.optim as optim
import yaml
from confgen.e2c.dataset import CustomData
from confgen.model.gnn import GNN
from confgen.molecule.graph import rdk2graph
from confgen.molecule.gt import isomorphic_core
from confgen.utils.utils import set_rdmol_positions, WarmCosine
from easydict import EasyDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.rdmolops import RemoveHs
from torch.optim.lr_scheduler import LambdaLR
from torch_geometric.data import DataLoader
from tqdm import tqdm


def featurize_mol_from_smiles(smiles, remove_hs=True):
    # filter fragments
    if '.' in smiles:
        return None

    # filter mols rdkit can't intrinsically handle
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        mol = Chem.AddHs(mol)
    else:
        return None
    N = mol.GetNumAtoms()

    conf = Chem.Conformer(N)
    for atom_id in range(N):
        conf.SetAtomPosition(atom_id, [0, 0, 0])
    mol.AddConformer(conf, assignId=True)

    if remove_hs:
        try:
            new_mol = RemoveHs(mol)
        except Exception:
            pass
    else:
        new_mol = mol

    graph = rdk2graph(new_mol)

    assert len(graph["edge_attr"]) == graph["edge_index"].shape[1]
    assert len(graph["node_feat"]) == graph["num_nodes"]

    data = CustomData()
    data.edge_index = torch.from_numpy(graph["edge_index"]).to(torch.int64)
    data.edge_attr = torch.from_numpy(graph["edge_attr"]).to(torch.int64)
    data.x = torch.from_numpy(graph["node_feat"]).to(torch.int64)
    data.n_nodes = graph["n_nodes"]
    data.n_edges = graph["n_edges"]

    data.rd_mol = copy.deepcopy(new_mol)
    data.isomorphisms = isomorphic_core(new_mol)

    data.nei_src_index = torch.from_numpy(graph["nei_src_index"]).to(torch.int64)
    data.nei_tgt_index = torch.from_numpy(graph["nei_tgt_index"]).to(torch.int64)
    data.nei_tgt_mask = torch.from_numpy(graph["nei_tgt_mask"]).to(torch.bool)

    data.pos = torch.zeros(data.n_nodes, 3, dtype=torch.float)

    return data


def evaluate_one(model, device, loader):
    model.eval()
    mol_preds = []
    for batch in tqdm(loader, desc="Iteration"):
        batch = batch.to(device)
        with torch.no_grad():
            pred, _ = model(batch)
        pred = pred[-1]
        batch_size = batch.num_graphs
        n_nodes = batch.n_nodes.tolist()
        pre_nodes = 0
        for i in range(batch_size):
            mol_preds.append(set_rdmol_positions(batch.rd_mol[i], pred[pre_nodes: pre_nodes + n_nodes[i]]))
            pre_nodes += n_nodes[i]

    mol = mol_preds[0]
    for m in mol_preds[1:]:
        mol.AddConformer(m.GetConformer(0), assignId=True)

    rdMolAlign.AlignMolConformers(mol)
    mol = Chem.AddHs(mol, addCoords=True)

    return mol


def repeat_data(data, num_repeat):
    data_list = [data.clone() for i in range(num_repeat)]
    return data_list


def dmcg_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    config_file = 'script/dmcg/config.yml'

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    np.random.seed(config.args.seed)
    torch.manual_seed(config.args.seed)
    torch.cuda.manual_seed(config.args.seed)
    random.seed(config.args.seed)

    shared_params = config.shared_params
    model = GNN(**shared_params).to(device)

    checkpoint = torch.load(config.args.eval_from, map_location=device)['model_state_dict']
    cur_state_dict = model.state_dict()
    del_keys = []

    for k in checkpoint.keys():
        if k not in cur_state_dict:
            del_keys.append(k)

    for k in del_keys:
        del checkpoint[k]

    model.load_state_dict(checkpoint)

    if config.args.use_adamw:
        dmcg_optimizer = optim.AdamW(model.parameters(), lr=config.args.lr, betas=(0.9, config.args.beta2), weight_decay=config.args.weight_decay)
    else:
        dmcg_optimizer = optim.Adam(model.parameters(), lr=config.args.lr, betas=(0.9, config.args.beta2), weight_decay=config.args.weight_decay)

    if not config.args.lr_warmup:
        scheduler = LambdaLR(dmcg_optimizer, lambda x: 1.0)
    else:
        lrscheduler = WarmCosine(tmax=len(train_loader) * config.args.period, warmup=int(4e3))
        scheduler = LambdaLR(dmcg_optimizer, lambda x: lrscheduler.step(x))

    input_lines = pd.read_csv(input_file, sep=',', skiprows=14, header=0).values
    mol_chunks = [input_lines[i:i + 100] for i in range(0, len(input_lines), 100)]

    writer = Chem.SDWriter(output_file_sdf)

    for chunk_index, chunk in enumerate(mol_chunks):
        pickle_slice = {}  # 类似这样 {'mol1': rdkit.Chem.rdchem.Mol, 'mol2': rdkit.Chem.rdchem.Mol}
        output_file_pkl = os.path.join(work_path, f'fastsmcg_output_slice{chunk_index:0>6}.pkl')
        for index, (smiles, title, num_conf) in enumerate(chunk):
            mol_title = f'{title}_{index + 1 + chunk_index * 100}'
            try:
                data = featurize_mol_from_smiles(smiles, remove_hs=True)
                data_list = repeat_data(data, num_conf)
                data_loader = DataLoader(dataset=data_list, batch_size=config.args.batch_size, shuffle=False, num_workers=config.args.num_workers)

                mol = evaluate_one(model, device, data_loader)
                mol.SetProp('_Name', mol_title)

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
    dmcg_gen(work_path, optimizer)
