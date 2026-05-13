#!/bin/bash
arboreto_with_multiprocessing.py \
so.renamed_res.0.2.loom \
mm_mgi_tfs.txt \
--method grnboost2 \
--output adj.tsv \
--num_workers 32 \
--seed 777