#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#===========================================================
# File Name: dmcg_gen.py
# Author: wangzhe
# E-mail: wangzhehyd@163.com
# Date: 2022-07-26 14:22:38
# Last Modified: 2023-01-28 15:17:11
# Version: 1.0
# Description: 
# Copyright: Hou group
#===========================================================

import os
import sys
import copy
import collections
import numpy as np
import random
import pickle
import pandas as pd
import argparse
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import RemoveHs
import torch
import torch.optim as optim
from torch.optim.lr_scheduler import LambdaLR
from torch_geometric.data import Data, DataLoader
from confgen.molecule.graph import rdk2graph
from confgen.molecule.gt import isomorphic_core
from confgen.e2c.dataset import CustomData 
from confgen.model.gnn import GNN
from confgen.utils.utils import set_rdmol_positions
from confgen.utils.utils import (
    WarmCosine,
    set_rdmol_positions,
    get_best_rmsd,
    evaluate_distance,
)


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
        conf.SetAtomPosition(atom_id, [0,0,0])
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
    mol_dict = {}
    for batch in tqdm(loader, desc="Iteration"):
        batch = batch.to(device)
        with torch.no_grad():
            pred, _ = model(batch)
        pred = pred[-1]
        batch_size = batch.num_graphs
        n_nodes = batch.n_nodes.tolist()
        pre_nodes = 0
        for i in range(batch_size):
            mol_preds.append(set_rdmol_positions(batch.rd_mol[i], pred[pre_nodes : pre_nodes + n_nodes[i]]))
            pre_nodes += n_nodes[i]
            #mol_dict[batch.title[i]] = mol_preds
    
    #return mol_dict
    return mol_preds

def repeat_data(data, num_repeat):
    datas = [data.clone() for i in range(num_repeat)]
    return datas

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--global-reducer", type=str, default="sum")
    parser.add_argument("--node-reducer", type=str, default="sum")
    parser.add_argument("--graph-pooling", type=str, default="sum")
    parser.add_argument("--dropedge-rate", type=float, default=0.1)
    parser.add_argument("--dropnode-rate", type=float, default=0.1)
    parser.add_argument("--num-layers", type=int, default=12)
    parser.add_argument("--decoder-layers", type=int, default=None)
    parser.add_argument("--latent-size", type=int, default=256)
    parser.add_argument("--mlp-hidden-size", type=int, default=1024)
    parser.add_argument("--mlp_layers", type=int, default=2)
    parser.add_argument("--use-layer-norm", action="store_true", default=False)
    parser.add_argument("--smiles", type=str, default="")
    parser.add_argument("--eval_from", type=str, default="")

    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--num-workers", type=int, default=1)
    # parser.add_argument("--log-dir", type=str, default="", help="tensorboard log directory")
    parser.add_argument("--checkpoint-dir", type=str, default="")
    parser.add_argument("--input_file", type=str, default="")
    parser.add_argument("--out", type=str, default="")

    parser.add_argument("--log-interval", type=int, default=100)
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--encoder-dropout", type=float, default=0.0)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--layernorm-before", action="store_true", default=False)
    parser.add_argument("--use-bn", action="store_true", default=False)
    parser.add_argument("--weight-decay", type=float, default=1e-2)
    parser.add_argument("--use-adamw", action="store_true", default=False)
    parser.add_argument("--beta2", type=float, default=0.999)
    parser.add_argument("--period", type=float, default=10)

    parser.add_argument("--base-path", type=str, default="")
    parser.add_argument("--dataset-name", type=str, default="qm9", choices=["qm9", "drugs", "iso17"])
    parser.add_argument("--train-size", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=2021)
    parser.add_argument("--lr-warmup", action="store_true", default=False)
    parser.add_argument("--enable-tb", action="store_true", default=False)
    parser.add_argument("--aux-loss", type=float, default=0.0)
    parser.add_argument("--train-subset", action="store_true", default=False)
    parser.add_argument("--eval-from", type=str, default=None)
    parser.add_argument("--data-split", type=str, choices=["cgcf", "default", "confgf"], default="default")
    parser.add_argument("--reuse-prior", action="store_true", default=False)
    parser.add_argument("--cycle", type=int, default=1)

    parser.add_argument("--vae-beta", type=float, default=1.0)
    parser.add_argument("--eval-one", action="store_true", default=False)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--extend-edge", action="store_true", default=False)
    parser.add_argument("--use-ff", action="store_true", default=False)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--pred-pos-residual", action="store_true", default=False)
    parser.add_argument("--node-attn", action="store_true", default=False)
    parser.add_argument("--global-attn", action="store_true", default=False)
    parser.add_argument("--shared-decoder", action="store_true", default=False)
    parser.add_argument("--shared-output", action="store_true", default=False)
    parser.add_argument("--sample-beta", type=float, default=1.0)
    parser.add_argument("--remove-hs", action="store_true", default=False)
    parser.add_argument("--use-global", action="store_true", default=False)
    parser.add_argument("--prop-pred", action="store_true", default=False)

    parser.add_argument("--score", action="store_true", default=False)
    parser.add_argument("--sigma-begin", type=float, default=10.0)
    parser.add_argument("--sigma-end", type=float, default=0.01)
    parser.add_argument("--noise-level", type=int, default=10)
    parser.add_argument("--noise-steps", type=int, default=100)
    parser.add_argument("--noise-lr", type=float, default=2.4e-6)
    parser.add_argument("--decoder-std", type=float, default=1.0)
    parser.add_argument("--score-prior", action="store_true", default=False)

    args = parser.parse_args()
    print(args)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed(args.seed)
    random.seed(args.seed)

    device = (
        torch.device("cuda:" + str(args.device))
        if torch.cuda.is_available()
        else torch.device("cpu")
    )

    shared_params = {
        "mlp_hidden_size": args.mlp_hidden_size,
        "mlp_layers": args.mlp_layers,
        "latent_size": args.latent_size,
        "use_layer_norm": args.use_layer_norm,
        "num_message_passing_steps": args.num_layers,
        "global_reducer": args.global_reducer,
        "node_reducer": args.node_reducer,
        "dropedge_rate": args.dropedge_rate,
        "dropnode_rate": args.dropnode_rate,
        "dropout": args.dropout,
        "layernorm_before": args.layernorm_before,
        "encoder_dropout": args.encoder_dropout,
        "use_bn": args.use_bn,
        "vae_beta": args.vae_beta,
        "decoder_layers": args.decoder_layers,
        "reuse_prior": args.reuse_prior,
        "cycle": args.cycle,
        "pred_pos_residual": args.pred_pos_residual,
        "node_attn": args.node_attn,
        "global_attn": args.global_attn,
        "shared_decoder": args.shared_decoder,
        "sample_beta": args.sample_beta,
        "shared_output": args.shared_output,
    }

    print(shared_params)
    model = GNN(**shared_params).to(device)
    if args.eval_from is not None:
        assert os.path.exists(args.eval_from)
        checkpoint = torch.load(args.eval_from, map_location=device)["model_state_dict"]
        print(len(checkpoint))
        cur_state_dict = model.state_dict()
        del_keys = []
        #print('cur_stat_dict: ', cur_state_dict)
        for k in checkpoint.keys():
            if k not in cur_state_dict:
                del_keys.append(k)
        print('del_keys: ', del_keys)
        for k in del_keys:
            del checkpoint[k]
        model.load_state_dict(checkpoint)

    num_params = sum(p.numel() for p in model.parameters())
    #print(f"#Params: {num_params}")

    if args.use_adamw:
        optimizer = optim.AdamW(model.parameters(), lr=args.lr, betas=(0.9, args.beta2), weight_decay=args.weight_decay)
    else:
        optimizer = optim.Adam(model.parameters(), lr=args.lr, betas=(0.9, args.beta2), weight_decay=args.weight_decay)

    if not args.lr_warmup:
        scheduler = LambdaLR(optimizer, lambda x: 1.0)
    else:
        lrscheduler = WarmCosine(tmax=len(train_loader) * args.period, warmup=int(4e3))
        scheduler = LambdaLR(optimizer, lambda x: lrscheduler.step(x))

    inputs = pd.read_csv(args.input_file)
    out = collections.OrderedDict()
    #for smiles, title, num_confs in tqdm(inputs.values):
    for smiles, num_confs, cs in tqdm(inputs.values):
        data = featurize_mol_from_smiles(cs, remove_hs=True)
        #data['title'] = title
        data['num_confs'] = num_confs
        data_list = repeat_data(data, num_confs*2) 

        data_loader = DataLoader(
            dataset=data_list,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=args.num_workers,
        )

        mol_preds = evaluate_one(model, device, data_loader)
        out[cs] = [ Chem.AddHs(mol, addCoords=True) for mol in mol_preds ]
 
    with open(f'{args.out}', 'wb') as f:
        pickle.dump(out, f)


if __name__ == '__main__':
    main()
