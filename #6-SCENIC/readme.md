# SCENIC使用

在服务器中，激活环境：

```bash
conda activate pyscenic
```

在服务器的SCENIC文件夹中，启动jupyter：

```bash
jupyter notebook --no-browser --port=8080
```

在本地powershell中，输入命令：

```bash
ssh -L 8080:localhost:8080 hekwlab@192.168.221.161
```

在本地浏览器中输入：

```bash
http://localhost:8080/
```

## Step 1: grn

GRN: Gene Regulatory Network, 基因调控网络

GRNBoost2: 一种基于梯度提升的有效调控网络推理算法，基于GENIE3，识别并筛选出与转录因子（TF）共表达的基因

GENIE3：Gene Network Inference with Ensemble of trees，使用基于树的集成学习算法

> 集成学习是指：组合多个相对较弱的学习算法以期获得更好的性能。
> 
> GENIE3内置了两种基于树的集成学习算法，Random Forests随机森林和Extra-Trees。
> 
> 随机森林：由多个决策树构成
> 
> 决策树举例：
> 
> ![10610675-37ad4eccde9ee464.png](C:\Users\admin\Desktop\10610675-37ad4eccde9ee464.png)
> 
> 基因的调控关系采用回归的方式获取。回归用于预测连续的、具体的数值。
> 
> ![10610675-0e0e0167d7976add.png](C:\Users\admin\Desktop\10610675-0e0e0167d7976add.png)
> 
> 上图中，gene1,2,3是转录因子， gene4是受gene1,2,3调控的靶基因。从上图决策树可以看出，gene4的表达受gene1,2,3的调控作用。它是一个很简单的回归树。gene1,2,3对gene4的调控作用不一，从而可以得到不同重要性分数。
> 
> 上图回归树思想为GENIE3用于构建调控网络关系的重要原因之一。当然作者没直接用回归决策树，而是使用基于决策树的集成学习算法，当然是因为决策树有缺点了。不过基本调控关系的构建思想还是回归树了。

大数据库用`arboreto_with_multiprocessing.py`脚本，sbatch到slurm服务器上：

在`pyscenic_arboreto.sh`文件中：

```bash
#!/bin/bash

# SBATCH -p short    # partition name
#SBATCH --nodelist=cpu04
#SBATCH -c 30                   # Request N cores [>12 cores recommended for pyscenic2 steps]
#SBATCH -t 5-00:00               # hours:minutes time after which job will be killed [16h needed for 50k cells x 20k cells with 20 cores]
#SBATCH --mem 150G              # total amount of memory requested [>100G recommended]
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written


# Move to desired working directory
cd /share/home/hekwlab/ZZY/SCENIC

# Set up Python environment
# module load anaconda3
  ## This module is necessary to use *source* activate pyscenic2
# conda activate pyscenic

# Double-check what Python version is running
python3 --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Create results directory
res_dir=outs_$(date +"%Y-%m-%d")
mkdir -p $res_dir
date

# Run pySCENIC get regulatory network from command line interface
# This circumvents the use of dask, which seems to cause trouble...
python3 ./arboreto_with_multiprocessing.py \
      /share/home/hekwlab/ZZY/SCENIC/so.renamed_res.0.2.loom \
        ./resources/allTFs_mm.txt \
        --num_workers 30 \
        -o $res_dir/pySCENIC_GRN_adjacencies.csv \
        --method grnboost2 \
        --sparse \
        --seed 787878

date
echo "regulatory networks inferred"
```

在linux中激活环境：

```bash
conda activate pyscenic
```

向slurm服务器提交任务：

```bash
sbatch pyscenic_arboreto.sh
```

## Step 2: Define modules of TF-regulons

run `pyscenic ctx`

This step cleans up the identified co-expression networks by drawing information from the feather databases and motif annotation files to identify target genes with significantly enriched motifs.

In `pyscenic_ctx.sh`:

```bash
#!/bin/bash

# SBATCH -p short    # partition name
#SBATCH --nodelist=cpu04
#SBATCH -c 30                   # Request N cores [>12 cores recommended for pyscenic2 steps]
#SBATCH -t 5-00:00               # hours:minutes time after which job will be killed [16h needed for 50k cells x 20k cells with 20 cores]
#SBATCH --mem 150G              # total amount of memory requested [>100G recommended]
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written


# Move to desired working directory
cd /share/home/hekwlab/ZZY/SCENIC

# Set up Python environment
# module load anaconda3
  ## This module is necessary to use *source* activate pyscenic2
# conda activate pyscenic

# Double-check what Python version is running
python3 --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH

# Clean up modules using feather ranking databases
pyscenic ctx outs_2023-04-12/pySCENIC_GRN_adjacencies.csv \
    databases/*.feather \
    --annotations_fname resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname /share/home/hekwlab/ZZY/SCENIC/so.renamed_res.0.2.loom \
    --output outs_2023-04-12/pySCENIC_CTX_regulons.csv \
    --mask_dropouts \
    --num_workers 30

date
echo "regulon modules defined"
sleep 5
```

## Step 3: Calculate regulon modules score for each cell

run `pyscenic aucell` to calculate a regulon module score for each cell in the dataset.

In `pyscenic_aucell.sh`:

```bash
#!/bin/bash

#SBATCH -c 12                   # Request N cores (must match future plan requested in R script)
#SBATCH -t 0-10:00               # hours:minutes runlimit after which job will be killed
#SBATCH --mem 100G               # total amount of memory requested
#SBATCH --job-name henderson    # Job name
#SBATCH -o %j.out               # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written

# Double-check what Python version is running
python --version
echo "PATH ="; echo $PATH
echo "PYTHONPATH ="; echo $PYTHONPATH


# Calculate cell module scores (from filtered loom, to extract acurate regulon incidence matrix later on)
pyscenic aucell \
    /share/home/hekwlab/ZZY/SCENIC/so.renamed_res.0.2.loom \
    ./outs_2023-04-12/pySCENIC_CTX_regulons.csv \
    --output ./outs_2023-04-12/pySCENIC_filtered.loom \
    --num_workers 12

echo "cell scoring completed"

# Add dimensionality reduction (derived from pySCENIC AUCell matrix)
   ## I always skipped this step, as it can also be performed in R (with proper sample integration)
python add_visualization.py \
   --loom_input ./outs_2023-04-12/pySCENIC_filtered.loom \
   --loom_output ./outs_2023-04-12/pySCENIC_viz.loom \
   --num_workers 12

#echo "dimensionality reduction added"
```
