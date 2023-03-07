import copy
import os
import pickle
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import torch
import yaml
from easydict import EasyDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.rdchem import HybridizationType, BondType
from torch_geometric.data import Batch
from torch_geometric.data import Data, Dataset
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
        # print('[Packed] %d Molecules, %d Conformations.' % (len(self._packed_data), len(self.data)))

        new_data = []
        # logic
        # save graph structure for each mol once, but store all confs
        cnt = 0
        for k, v in self._packed_data.items():
            data = copy.deepcopy(v[0])
            all_pos = []
            for i in range(len(v)):
                all_pos.append(v[i].pos)
            data.pos_ref = torch.cat(all_pos, 0)  # (num_conf*num_node, 3)
            data.num_pos_ref = torch.tensor([len(all_pos)], dtype=torch.long)
            # del data.pos

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


def repeat_data(data, num_repeat):
    datas = [data.clone() for i in range(num_repeat)]
    return Batch.from_data_list(datas)


def geodiff_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    config_file = 'script/geodiff/config.yml'
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)

    gpus = list(filter(lambda x: x is not None, config.train.gpus))
    device = torch.device(gpus[0]) if len(gpus) > 0 else torch.device('cpu')
    config.train.device = device
    config.train.gpus = gpus

    seed = config.train.seed
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)


    ckpt = torch.load(f'{config.test.init_checkpoint}/{config.test.epoch}.pt')
    model = get_model(ckpt['config'].model).to(config.train.device)
    model.load_state_dict(ckpt['model'])

    transforms = Compose([CountNodesPerGraph(), AddHigherOrderEdges(order=config.model.edge_order)])

    input_lines = pd.read_csv(input_file, sep=',', skiprows=14, header=0).values
    mol_chunks = [input_lines[i:i + 100] for i in range(0, len(input_lines), 100)]

    writer = Chem.SDWriter(output_file_sdf)
    for chunk_index, chunk in enumerate(mol_chunks):
        pickle_slice = {}  # 类似这样 {'mol1': rdkit.Chem.rdchem.Mol, 'mol2': rdkit.Chem.rdchem.Mol}
        output_file_pkl = os.path.join(work_path, f'fastsmcg_output_slice{chunk_index:0>6}.pkl')
        for index, (smiles, title, num_conf) in enumerate(chunk):
            mol_title = f'{title}_{index + 1 + chunk_index * 100}'
            try:
                data_init = smiles_to_data(smiles)
                data_init['title'] = mol_title
                data_init['num_confs'] = num_conf
                dataset = PackedConformationDataset([data_init], transform=transforms)

                data = dataset[0]
                data_input = data.clone()
                data_input['pos_ref'] = None
                batch = repeat_data(data_input, num_conf).to(config.train.device)
                clip_local = None
                for _ in range(2):  # Maximum number of retry
                    try:
                        pos_init = torch.randn(batch.num_nodes, 3).to(config.train.device)
                        pos_gen, pos_gen_traj = model.langevin_dynamics_sample(
                            atom_type=batch.atom_type,
                            pos_init=pos_init,
                            bond_index=batch.edge_index,
                            bond_type=batch.edge_type,
                            batch=batch.batch,
                            num_graphs=batch.num_graphs,
                            extend_order=False,  # Done in transforms.
                            n_steps=5000,
                            step_lr=1e-6,
                            w_global=1.0,
                            global_start_sigma=0.5,
                            clip=1000.0,
                            clip_local=clip_local,
                            sampling_type='ld',
                            eta=1.0
                        )
                        pos_gen = pos_gen.cpu()

                        if config.test.save_traj:
                            data.pos_gen = torch.stack(pos_gen_traj)
                        else:
                            data.pos_gen = pos_gen

                        # 返回的数据处理成 rdkit.Chem.rdchem.Mol
                        mol = data.rdmol
                        mol.SetProp('_Name', mol_title)
                        num_confs = int(data.pos_gen.shape[0] / data.pos.shape[0])
                        pos_gen_list = data.pos_gen.tolist()
                        coords_list = [pos_gen_list[i:i + data.pos.shape[0]] for i in range(0, data.pos_gen.shape[0], data.pos.shape[0])]

                        for conf_id in range(num_confs):
                            conf = Chem.Conformer(mol.GetNumAtoms())
                            for atom_id in range(mol.GetNumAtoms()):
                                conf.SetAtomPosition(atom_id, coords_list[conf_id][atom_id])
                            mol.AddConformer(conf, assignId=True)
                        rdMolAlign.AlignMolConformers(mol)
                        break  # No errors occured, break the retry loop
                    except:
                        clip_local = 20

                if mol.GetNumConformers() > 0:
                    # 优化构象
                    for conf_id in range(mol.GetNumConformers()):
                        if optimizer == '1':
                            AllChem.UFFOptimizeMolecule(mol, maxIters=500, confId=conf_id)
                        elif optimizer == '2':
                            AllChem.MMFFOptimizeMolecule(mol, maxIters=500, confId=conf_id)
                        elif optimizer == '3':
                            AllChem.MMFFOptimizeMolecule(mol, maxIters=500, confId=conf_id, mmffVariant='MMFF94s')
                        writer.write(mol, confId=conf_id)
                        pickle_slice[mol_title] = mol

                else:
                    pickle_slice[mol_title] = []
            except:
                pickle_slice[mol_title] = []

        with open(output_file_pkl, 'wb') as f:
            pickle.dump(pickle_slice, f)

    writer.flush()
    writer.close()



if __name__ == '__main__':
    work_path = sys.argv[1]
    optimizer = sys.argv[2]
    geodiff_gen(work_path, optimizer)
