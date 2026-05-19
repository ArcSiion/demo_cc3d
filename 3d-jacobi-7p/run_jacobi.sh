#!/bin/bash
# ==============================
# BitLab 集群演示脚本：run_jacobi.sh
# 用途：演示如何用 sbatch 提交一个纯 CPU 的批处理任务
# 说明：本脚本不申请 GPU，不激活 conda，仅执行当前目录下的 test.sh
# ==============================


# ==================================================
# 区域一：资源申请（告诉 Slurm 你需要多少算力）
# ==================================================

# 作业名称：提交后在 squeue 里能看到这个名字
#SBATCH --job-name=jacobi_demo

# [备份] 原 GPU 申请方式，演示时不需要 GPU，先注释保留
##SBATCH --gres=gpu:1

# 申请 4 个 CPU 核心
#SBATCH --cpus-per-task=4

# 申请 8G 内存
#SBATCH --mem=8G

# 最长运行时间设置为 10 分钟
#SBATCH --time=00:10:00

# 邮件通知：任务开始、结束、失败时都发送邮件
#SBATCH --mail-type=ALL
#SBATCH --mail-user=18551224753@163.com

# ==================================================
# 以下参数为实验室统一格式，通常不需要改
# ==================================================

# 指定队列
#SBATCH --partition=compute

# 强制单节点运行
#SBATCH --nodes=1

# 单任务模式
#SBATCH --ntasks=1

# 标准输出日志文件
#SBATCH --output=log/run_%j.log

# 标准错误日志文件
#SBATCH --error=log/err_%j.log


# ==================================================
# 区域二：运行环境（说明任务在哪个目录执行）
# ==================================================

# 切换到提交作业时所在的目录
# 这样无论你从哪个位置 sbatch，本脚本都会回到提交目录执行
cd "$SLURM_SUBMIT_DIR"

echo "=== BitLab 算力集群: 演示任务开始 ==="
echo "作业ID: $SLURM_JOB_ID"
echo "节点名: $SLURMD_NODENAME"
echo "工作目录: $SLURM_SUBMIT_DIR"
echo "分配的 CPU 核心数: $SLURM_CPUS_PER_TASK"

# 将 OpenMP 线程数设置为和申请到的 CPU 核心数一致
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# [备份] 原 conda 环境激活方式，演示时不需要，先注释保留
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate my_pytorch_env


# ==================================================
# 区域三：执行代码（真正运行你的程序）
# ==================================================

echo "=== 开始执行 test.sh ==="
bash test.sh
echo "=== BitLab 算力集群: 演示任务结束 ==="
