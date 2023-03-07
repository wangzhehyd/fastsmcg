import copy
import os
import pickle
import random
import sys

import numpy as np
import pandas as pd
import torch
import yaml
from confgf import models, dataset, runner, utils
from easydict import EasyDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


# 继承runner中的DefaultRunner类 重写generate_samples_from_smiles方法
class SMCGRunner(runner.DefaultRunner):
    def generate_samples_from_smiles(self, smiles, generator, num_repeat=1, keep_traj=False):
        data = dataset.smiles_to_data(smiles)

        if data is None:
            raise ValueError('invalid smiles: %s' % smiles)

        return_data = copy.deepcopy(data)
        batch = utils.repeat_data(data, num_repeat).to(self.device)

        if generator == 'ConfGF':
            _, _, batch, pos_traj = self.ConfGF_generator(batch, self.config.test.gen)  # (sigams, 100, num_node, 3)
        elif generator == 'ConfGFDist':
            embedder = utils.Embed3D(step_size=self.config.test.gen.dg_step_size, num_steps=self.config.test.gen.dg_num_steps, verbose=self.config.test.gen.verbose)
            _, _, batch, pos_traj = self.ConfGFDist_generator(batch, self.config.test.gen, embedder)  # (num_steps, num_node, 3)
        else:
            raise NotImplementedError

        batch = batch.to('cpu').to_data_list()
        pos_traj = pos_traj.view(-1, 3).to('cpu')
        pos_traj_step = pos_traj.size(0) // return_data.num_nodes

        all_pos = []
        for i in range(len(batch)):
            all_pos.append(batch[i].pos_gen)
        return_data.pos_gen = torch.cat(all_pos, 0)  # (num_repeat * num_node, 3)
        return_data.num_pos_gen = torch.tensor([len(all_pos)], dtype=torch.long)
        if keep_traj:
            return_data.pos_traj = pos_traj
            return_data.num_pos_traj = torch.tensor([pos_traj_step], dtype=torch.long)

        # 返回的数据处理成 rdkit.Chem.rdchem.Mol
        mol = return_data.rdmol
        num_confs = int(return_data.pos_gen.shape[0] / return_data.pos.shape[0])

        pos_gen_list = return_data.pos_gen.tolist()
        coords_list = [pos_gen_list[i:i + return_data.pos.shape[0]] for i in range(0, return_data.pos_gen.shape[0], return_data.pos.shape[0])]

        for conf_id in range(num_confs):
            conf = Chem.Conformer(mol.GetNumAtoms())
            for atom_id in range(mol.GetNumAtoms()):
                conf.SetAtomPosition(atom_id, coords_list[conf_id][atom_id])
            mol.AddConformer(conf, assignId=True)
        rdMolAlign.AlignMolConformers(mol)

        return mol


def confgf_gen(work_path, optimizer):
    input_file = os.path.join(work_path, 'fastsmcg_input.txt')
    output_file_sdf = os.path.join(work_path, 'fastsmcg_output.sdf')

    config_file = 'script/confgf/config.yml'

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config = EasyDict(config)

    gpus = list(filter(lambda x: x is not None, config.train.gpus))
    device = torch.device(gpus[0]) if len(gpus) > 0 else torch.device('cpu')
    config.train.device = device
    config.train.gpus = gpus

    seed = config.train.seed
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = True

    model = models.DistanceScoreMatch(config)
    solver = SMCGRunner(train_set=None, val_set=None, test_set=None, model=model, optimizer=None, scheduler=None, gpus=gpus, config=config)
    solver.load(config.test.init_checkpoint, epoch=config.test.epoch)

    input_lines = pd.read_csv(input_file, sep=',', skiprows=14, header=0).values
    mol_chunks = [input_lines[i:i + 100] for i in range(0, len(input_lines), 100)]

    writer = Chem.SDWriter(output_file_sdf)

    for chunk_index, chunk in enumerate(mol_chunks):
        pickle_slice = {}  # 类似这样 {'mol1': rdkit.Chem.rdchem.Mol, 'mol2': rdkit.Chem.rdchem.Mol}
        output_file_pkl = os.path.join(work_path, f'fastsmcg_output_slice{chunk_index:0>6}.pkl')
        for index, (smiles, title, num_conf) in enumerate(chunk):
            mol_title = f'{title}_{index + 1 + chunk_index * 100}'
            try:
                mol = solver.generate_samples_from_smiles(smiles, 'ConfGF', num_repeat=num_conf, keep_traj=False)
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
    confgf_gen(work_path, optimizer)
