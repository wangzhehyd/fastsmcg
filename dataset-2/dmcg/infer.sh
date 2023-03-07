. ~/.conda/setenv.sh && conda activate dmcg && python dmcg_gen.py --dropout 0.1 \
    --input_file ../dataset-2.csv --use-bn  --use-adamw --train-subset --num-layers 6 \
    --eval-from checkpoint_94.pt --workers 20 --batch-size 128 \
    --reuse-prior --data-split default --node-attn --dataset-name drugs \
    --remove-hs --shared-output --pred-pos-residual --sample-beta 1.2 \
    --mlp-hidden-size 1024 --latent-size 256 --num-workers 0 --out dataset-2-out.pkl
