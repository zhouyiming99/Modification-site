#!/bin/bash

# --------------------
# 【依赖版本要求】
# java openjdk    1.8.0
# minimap2        2.14-r886
# samtools        0.1.19
# sam2tsv         a779a30d6af485d9cd669aa3752465132cf21eec
# python          3.6.7
# h5py            2.8.0
# numpy           1.15.4
# pandas          0.23.4
# scikit-learn    0.20.2
# nanopolish      0.12.4
# dask            2.5.2
# biopython       1.76
# pysam           0.15.3+
# R               3.6.0
# R packages: forcats / optparse / stringr / dplyr / purrr /
#             readr / tidyr / tibble / tidyverse / ggrepel /
#             car / ggplot2 / reshape2 / outliers
# --------------------

# --------------------
# 【用户参数配置区】
# --------------------

# EpiNano 软件安装目录
Epinano=path_to_Epinano_software

# 参考转录组 FASTA 文件
ref=path_to_reference_transcriptome

# 经过 minimap2 比对并排序后的 BAM 文件
bam=path_to_bam_file

# 样本名前缀（用于输出文件命名）
SAMPLE="sample"

# --------------------
# 【前置步骤】minimap2 比对（如已有 BAM 文件可跳过）
# 将 fastq reads 比对到参考转录组，生成排序后的 BAM 文件
# -ax map-ont：Nanopore 长读长模式
# -uf：正链比对（直接 RNA 测序）
# --secondary=no：不输出次优比对结果
# --------------------

minimap2 \
    -ax map-ont \
    -uf \
    --secondary=no \
    -t 16 \
    ${ref} \
    ${SAMPLE}.fastq.gz \
    | samtools view -Sb \
    | samtools sort -o ${SAMPLE}.bam

# 对 BAM 文件建立索引
samtools index ${SAMPLE}.bam

# --------------------
# 【步骤 1】提取碱基变异特征
# 从 BAM 文件中逐位点统计碱基质量、错配率、缺失率、插入率等特征
# -R：参考转录组 FASTA 文件
# -b：输入 BAM 文件
# -s：sam2tsv 工具路径（Java 程序，用于解析 CIGAR 字段）
# -n：并行 CPU 线程数
# -T：输出类型（t 表示转录组模式）
# 输出文件：${SAMPLE}.plus_strand.per.site.csv
# --------------------
python $Epinano/Epinano_Variants.py \
    -R $ref \
    -b $bam \
    -s $Epinano/misc/sam2tsv.jar \
    -n 16 \
    -T t

# --------------------
# 【步骤 2】滑动窗口合并特征
# 将逐位点特征合并为 5-mer 窗口特征（以中心位点为基准，左右各取 2 个位点）
# 输出文件：${SAMPLE}.plus_strand.per.site.5mer.csv
# --------------------
python $Epinano/misc/Slide_Variants.py \
    ${SAMPLE}.plus_strand.per.site.csv \
    5                               # 滑动窗口大小（推荐使用 5-mer）

# --------------------
# 【步骤 3】预测 RNA 修饰位点
# 使用预训练的线性 SVM 模型对每个 5-mer 位点进行修饰概率预测
# --model：预训练模型路径（RRACH 基序 m6A 检测模型）
# --predict：输入 5-mer 特征文件
# --columns：用于预测的特征列索引（8=错配率, 13=插入率, 23=缺失率）
# --out_prefix：输出文件前缀
# 输出文件：${SAMPLE}_predictions.csv
# --------------------
python $Epinano/Epinano_Predict.py \
    --model $Epinano/models/rrach.q3.mis3.del3.linear.dump \
    --predict ${SAMPLE}.plus_strand.per.site.5mer.csv \
    --columns 8,13,23 \             # 质量分数均值、错配率、缺失率对应列索引
    --out_prefix ${SAMPLE}
