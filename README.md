# fastSMCG
Data and code required to reach the main conclusions of the benchmark paper (Small-Molecule Conformer Generators: Evaluation of Traditional Methods and AI Models on High-quality Datasets)

fastSMCG is a user-friendly webserver for small molecule conformer generation and analysis that integrates RDKit ETKDG algorithm and multiple AI models. It can be accessed free of charge at http://cadd.zju.edu.cn/fastsmcg/.

## dataset-1
### analysis directory
Contains the result raw data and demo scripts for plotting figures  

### script directory
Contains demo scripts for generating conformaers. For example, the ConfGenx method can be called to generate conformaers with a maximum ensemble size of 1 by running the command below  
```bash
sbatch confgenx_001_submit.slurm
```

### Raw structures of Dataset-I for benchmark
csv format: dataset-1.csv  
pickle format: dataset-1.pkl  
SD format: dataset.sdf  
SMILES format: dataset.smi  

### Raw structures of Dataset-I for benchmark （molecule title with index suffix）
csv format: dataset-1_with-index.csv  
pickle format: dataset-1_with-index.pkl  
SD format: dataset-1_with-index.sdf  
SMILES format: dataset-1_with-index.smi  
SD format: dataset-1_with-index_2d.sdf (contains 2D coordinates as input for ConfGenx method)  

### run_fastsmcg.sh
A demo script for calling AI models or RDKit ETKDG method to generate conformers  
the RDKit ETKDG method can be called to generate conformaers with a maximum ensemble size of 1 by running the command below 
```bash
bash run_fastsmcg.sh rdkit 1
```

## dataset-2
### analysis
Contains the result raw data and demo script for plotting figures  

### script directory
Contains demo scripts for generating conformaers. For example, the ConfGenx method can be called to generate conformaers by running the command below  
```bash
sbatch confgenx_submit.slurm
```

### Raw structures of Dataset-II for benchmark
csv format: dataset-2.csv  
pickle format: dataset-2.pkl 
