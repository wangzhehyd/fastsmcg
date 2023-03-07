#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: geodiff_gen.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-07-22 17:15:36
# Last Modified: 2022-10-28 16:43:53
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import copy
import argparse
import pickle
import random
import yaml
import numpy as np
import pandas as pd
from glob import glob
from collections import defaultdict
from tqdm.auto import tqdm
from easydict import EasyDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import HybridizationType, BondType
import torch
from torch_geometric.data import Data, Dataset
from torch_geometric.data import Batch
from torch_geometric.transforms import Compose
from models.epsnet import get_model
from utils.transforms import AddHigherOrderEdges

BOND_TYPES = {t: i for i, t in enumerate(BondType.names.values())}

def smiles_to_data(smiles):
    if '.' in smiles:
        return None
    else:
        try:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        except:
            return None

    N = mol.GetNumAtoms()
    pos = torch.rand((N, 3), dtype=torch.float32)

    atomic_number = []
    aromatic = []
    sp = []
    sp2 = []
    sp3 = []

    for atom in mol.GetAtoms():
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization = atom.GetHybridization()
        sp.append(1 if hybridization == HybridizationType.SP else 0)
        sp2.append(1 if hybridization == HybridizationType.SP2 else 0)
        sp3.append(1 if hybridization == HybridizationType.SP3 else 0)

    z = torch.tensor(atomic_number, dtype=torch.long)

    row, col, edge_type = [], [], []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [BOND_TYPES[bond.GetBondType()]]

    edge_index = torch.tensor([row, col], dtype=torch.long)
    edge_type = torch.tensor(edge_type)

    perm = (edge_index[0] * N + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]
    row, col = edge_index

    data = Data(atom_type=z, pos=pos, edge_index=edge_index, edge_type=edge_type, rdmol=copy.deepcopy(mol), smiles=smiles)

    return data

class ConformationDataset(Dataset):
    def __init__(self, data, transform=None):
        super().__init__()
        self.data = data
        self.transform = transform
        self.atom_types = self._atom_types()
        self.edge_types = self._edge_types()

    def __getitem__(self, idx):

        data = self.data[idx].clone()
        if self.transform is not None:
            data = self.transform(data)
        return data

    def __len__(self):
        return len(self.data)

    def _atom_types(self):
        """All atom types."""
        atom_types = set()
        for graph in self.data:
            atom_types.update(graph.atom_type.tolist())
        return sorted(atom_types)

    def _edge_types(self):
        """All edge types."""
        edge_types = set()
        for graph in self.data:
            edge_types.update(graph.edge_type.tolist())
        return sorted(edge_types)

class PackedConformationDataset(ConformationDataset):
    def __init__(self, path, transform=None):
        super().__init__(path, transform)
        self._pack_data_by_mol()

    def _pack_data_by_mol(self):
        """
        pack confs with same mol into a single data object
        """
        self._packed_data = defaultdict(list)
        if hasattr(self.data, 'idx'):
            for i in range(len(self.data)):
                self._packed_data[self.data[i].idx.item()].append(self.data[i])
        else:
            for i in range(len(self.data)):
                self._packed_data[self.data[i].smiles].append(self.data[i])
        print('[Packed] %d Molecules, %d Conformations.' % (len(self._packed_data), len(self.data)))

        new_data = []
        # logic
        # save graph structure for each mol once, but store all confs 
        cnt = 0
        for k, v in self._packed_data.items():
            data = copy.deepcopy(v[0])
            all_pos = []
            for i in range(len(v)):
                all_pos.append(v[i].pos)
            data.pos_ref = torch.cat(all_pos, 0) # (num_conf*num_node, 3)
            data.num_pos_ref = torch.tensor([len(all_pos)], dtype=torch.long)
            #del data.pos

            if hasattr(data, 'totalenergy'):
                del data.totalenergy
            if hasattr(data, 'boltzmannweight'):
                del data.boltzmannweight
            new_data.append(data)
        self.new_data = new_data

    def __getitem__(self, idx):

        data = self.new_data[idx].clone()
        if self.transform is not None:
            data = self.transform(data)
        return data

    def __len__(self):
        return len(self.new_data)

class CountNodesPerGraph(object):
    def __init__(self) -> None:
        super().__init__()

    def __call__(self, data):
        data.num_nodes_per_graph = torch.LongTensor([data.num_nodes])
        return data

def seed_all(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

def repeat_data(data, num_repeat):
    datas = [data.clone() for i in range(num_repeat)]
    return Batch.from_data_list(datas)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ckpt', type=str, help='path for loading the checkpoint')
    parser.add_argument('--save_traj', action='store_true', default=False, help='whether store the whole trajectory for sampling')
    parser.add_argument('--input_file', type=str, default=None)
    parser.add_argument('--out', type=str, default=None)
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--clip', type=float, default=1000.0)
    parser.add_argument('--n_steps', type=int, default=5000, help='sampling num steps; for DSM framework, this means num steps for each noise scale')
    parser.add_argument('--global_start_sigma', type=float, default=0.5, help='enable global gradients only when noise is low')
    parser.add_argument('--w_global', type=float, default=1.0, help='weight for global gradients')
    # Parameters for DDPM
    parser.add_argument('--sampling_type', type=str, default='ld', help='generalized, ddpm_noisy, ld: sampling method for DDIM, DDPM or Langevin Dynamics')
    parser.add_argument('--eta', type=float, default=1.0, help='weight for DDIM and DDPM: 0->DDIM, 1->DDPM')
    args = parser.parse_args()

    # Load checkpoint
    ckpt = torch.load(args.ckpt)
    config_path = glob(os.path.join(os.path.dirname(os.path.dirname(args.ckpt)), '*.yml'))[0]
    with open(config_path, 'r') as f:
        config = EasyDict(yaml.safe_load(f))
    seed_all(config.train.seed)

    # Datasets and loaders
    print('Loading datasets...')
    transforms = Compose([
        CountNodesPerGraph(),
        AddHigherOrderEdges(order=config.model.edge_order), # Offline edge augmentation
    ])

    inputs = pd.read_csv(args.input_file)
    all_data = []
    #for smiles, title, num_confs in tqdm(inputs.values):
    mol_index=1
    for smiles, title, num_confs in tqdm(inputs.values):
        data = smiles_to_data(smiles)
        #data['title'] = title
        data['title'] = f'mol_{mol_index}'
        data['num_confs'] = num_confs
        all_data.append(data)
        mol_index += 1

    dataset = PackedConformationDataset(all_data, transform=transforms)

    # Model
    print('Loading model...')
    model = get_model(ckpt['config'].model).to(args.device)
    model.load_state_dict(ckpt['model'])
    results = []

    for i, data in enumerate(tqdm(dataset)):
        if data.num_confs < 1000:
            num_samples = data.num_confs * 2
        
        if os.path.exists(f'geodiff_{data.title}_conf_{data.num_confs}.pkl'):
            continue

        data_input = data.clone()
        data_input['pos_ref'] = None
        batch = repeat_data(data_input, num_samples).to(args.device)
        
        clip_local = None
        for _ in range(2):  # Maximum number of retry
            try:
                pos_init = torch.randn(batch.num_nodes, 3).to(args.device)
                pos_gen, pos_gen_traj = model.langevin_dynamics_sample(
                    atom_type=batch.atom_type,
                    pos_init=pos_init,
                    bond_index=batch.edge_index,
                    bond_type=batch.edge_type,
                    batch=batch.batch,
                    num_graphs=batch.num_graphs,
                    extend_order=False, # Done in transforms.
                    n_steps=args.n_steps,
                    step_lr=1e-6,
                    w_global=args.w_global,
                    global_start_sigma=args.global_start_sigma,
                    clip=args.clip,
                    clip_local=clip_local,
                    sampling_type=args.sampling_type,
                    eta=args.eta
                )
                pos_gen = pos_gen.cpu()
                if args.save_traj:
                    data.pos_gen = torch.stack(pos_gen_traj)
                else:
                    data.pos_gen = pos_gen
                results.append(data)

                save_path = f'geodiff_{data.title}_conf_{data.num_confs}.pkl'
                with open(save_path, 'wb') as f:
                    pickle.dump([data], f)

                break   # No errors occured, break the retry loop
            except FloatingPointError:
                clip_local = 20
    save_path= f'{args.out}'
    with open(save_path, 'wb') as f:
        pickle.dump(results, f)
