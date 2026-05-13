#!/bin/bash
# sbatch -J test
# #指定任务队列，默认使用normal
# sbatch -p normal
# #申请4个节点，单节点任务可以省略不写
# sbatch -N 4
# #每个节点分配20个MPI进程，即cpu核心数。
# sbatch --ntasks-per-node=20
# #一共跑80个进程。也可以省略“-N”和” --ntasks-per-node”参数，这样任务会优先填满一个节点，不会平均分散在几个节点上。
#SBATCH --nodelist=cpu06
#SBATCH -n 100
#每个节点申请一个gpu，cpu任务请注释掉
# SBATCH --gres=gpu:1
# SBATCH --get-user-env
# ###脚本错误输出
# SBATCH -e job-%j.err
# ###脚本输出
# SBATCH -o job-%j.out
# ###进入提交任务时所在目录
# # cd $SLURM_SUBMIT_DIR
# ###加载环境变量，根据实际使用的程序修改
# module load openmpi
# module load miniconda3
####程序命令
Rscript SCENIC.R
date -R
# mpirun hostname
sleep 10
echo Completed!