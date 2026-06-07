#!/bin/bash


# bedtools      v2.29.2
# samtools      1.3.1
# minimap2      2.17-r941
# python        3.7.3
# h5py          2.9.0
# statsmodels   0.10.0
# joblib        0.16.0
# xgboost       0.80
# pysam         0.16.0.1
# tqdm          4.39.0
# pycairo       1.19.1
# scikit-learn  0.22
# --------------------


# Nanom6A 软件安装目录
Nanom6A=path_to_Nanom6A_software

# 经过 Tombo resquiggle 处理后的 fast5 文件目录
FAST5=path_to_basecalled_fast5_files

# 参考基因组 FASTA 文件
genome=path_to_reference_genome

# 参考基因组 BED 文件
# BED 格式如下（每列依次为：染色体、起始位置、终止位置、基因名、碱基、链方向）：
# chr1  3073252  3073253  RP23-271O17.1  A  +
bed=path_to_reference_genome_gene_bed

# --------------------
# 【步骤 1】列出所有 fast5 文件路径
# 递归搜索 workspace 目录下所有经过 Tombo resquiggle 的 fast5 文件
# --------------------
find $FAST5/workspace -name "*.fast5" > fast5.txt

# --------------------
# 【步骤 2】提取原始信号特征
# 从 fast5 文件中提取用于 m6A 预测的原始电信号特征
# --basecall_group：指定 Tombo resquiggle 写入的 HDF5 数据组名称
# --cpu：并行 CPU 线程数
# --clip：信号裁剪长度（0 表示不裁剪）
# --fl：输入 fast5 文件列表
# --------------------
mkdir -p results

python $Nanom6A/extract_raw_and_feature_fast.py \
    -o results/sample \
    --basecall_group RawGenomeCorrected_000 \
    --cpu 16 \
    --clip 0 \
    --fl fast5.txt

# --------------------
# 【步骤 3】预测 m6A 修饰位点
# 基于 XGBoost 模型对每个位点进行 m6A 修饰概率预测
# --support：每个位点最少支持的 reads 数
# --model：预训练模型路径
# --------------------
mkdir -p results_final

python $Nanom6A/predict_sites.py \
    -i results/sample \
    -o results_final/sample \
    -g $genome \
    -r $bed \
    --cpu 16 \
    --support 5 \
    --model $Nanom6A/bin/model
