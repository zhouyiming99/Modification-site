#!/bin/bash



# CHEUI 安装目录
CHEUI_PATH=~/CHEUI

# k-mer 参考模型文件（CHEUI 自带）
KMER_MODEL=${CHEUI_PATH}/kmer_models/model_kmer.csv

# 预训练深度学习模型目录（CHEUI 自带）
MODEL_DIR=${CHEUI_PATH}/CHEUI_trained_models

# nanopolish eventalign 输出文件（两个条件）
CONDITION_X_INPUT="nanopolish_X_eventalign.txt"    # 条件 X（对照组）
CONDITION_Y_INPUT="nanopolish_Y_eventalign.txt"    # 条件 Y（实验组）

# 每个位点最少覆盖 reads 数
MIN_READS=20

# 并行线程数
THREADS=8

# ====================
# ====== m6A 检测 ======
# ====================

# 从 eventalign 文件中提取以 A 为中心的 9-mer 信号（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_preprocess_m6A.py \
    -i ${CONDITION_X_INPUT} \
    -m ${KMER_MODEL} \
    -o condition_X_m6A_signals.p \
    -n ${THREADS}

# 从 eventalign 文件中提取以 A 为中心的 9-mer 信号（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_preprocess_m6A.py \
    -i ${CONDITION_Y_INPUT} \
    -m ${KMER_MODEL} \
    -o condition_Y_m6A_signals.p \
    -n ${THREADS}

# 模型 1 预测：计算每条 read 上每个 A 位点的 m6A 修饰概率（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model1.py \
    -i condition_X_m6A_signals.p \
    -m ${MODEL_DIR}/CHEUI_m6A_model1.h5 \
    -l X \
    -o condition_X_m6A_read_level.txt

# 模型 1 预测：计算每条 read 上每个 A 位点的 m6A 修饰概率（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model1.py \
    -i condition_Y_m6A_signals.p \
    -m ${MODEL_DIR}/CHEUI_m6A_model1.h5 \
    -l Y \
    -o condition_Y_m6A_read_level.txt

# 按转录本坐标排序（模型 2 要求输入文件有序）
sort -k1,1 -k2,2n condition_X_m6A_read_level.txt > condition_X_m6A_read_level_sorted.txt
sort -k1,1 -k2,2n condition_Y_m6A_read_level.txt > condition_Y_m6A_read_level_sorted.txt

# 模型 2 预测：汇总位点水平的 m6A 修饰概率与化学计量比（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model2.py \
    -i condition_X_m6A_read_level_sorted.txt \
    -m ${MODEL_DIR}/CHEUI_m6A_model2.h5 \
    -c 0 \                      # 输出所有位点（不设概率阈值过滤）
    -d 0.3 0.7 \                # read 分类双阈值（<0.3 未修饰，>0.7 修饰）
    -n ${MIN_READS} \
    -o condition_X_m6A_site_level.txt

# 模型 2 预测：汇总位点水平的 m6A 修饰概率与化学计量比（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model2.py \
    -i condition_Y_m6A_read_level_sorted.txt \
    -m ${MODEL_DIR}/CHEUI_m6A_model2.h5 \
    -c 0 \
    -d 0.3 0.7 \
    -n ${MIN_READS} \
    -o condition_Y_m6A_site_level.txt

# 生成差异分析配置文件，并使用 Mann-Whitney U 检验进行差异 m6A 分析
echo -e "X\tcondition_X_m6A_read_level_sorted.txt\nY\tcondition_Y_m6A_read_level_sorted.txt" \
    > m6A_diff_config.txt

python3 ${CHEUI_PATH}/scripts/CHEUI_diff.py \
    -i m6A_diff_config.txt \
    -n ${MIN_READS} \
    -o m6A_differential_results.txt

# ====================
# ====== m5C 检测 ======
# ====================

# 从 eventalign 文件中提取以 C 为中心的 9-mer 信号（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_preprocess_m5C.py \
    -i ${CONDITION_X_INPUT} \
    -m ${KMER_MODEL} \
    -o condition_X_m5C_signals.p \
    -n ${THREADS}

# 从 eventalign 文件中提取以 C 为中心的 9-mer 信号（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_preprocess_m5C.py \
    -i ${CONDITION_Y_INPUT} \
    -m ${KMER_MODEL} \
    -o condition_Y_m5C_signals.p \
    -n ${THREADS}

# 模型 1 预测：计算每条 read 上每个 C 位点的 m5C 修饰概率（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model1.py \
    -i condition_X_m5C_signals.p \
    -m ${MODEL_DIR}/CHEUI_m5C_model1.h5 \
    -l X \
    -o condition_X_m5C_read_level.txt

# 模型 1 预测：计算每条 read 上每个 C 位点的 m5C 修饰概率（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model1.py \
    -i condition_Y_m5C_signals.p \
    -m ${MODEL_DIR}/CHEUI_m5C_model1.h5 \
    -l Y \
    -o condition_Y_m5C_read_level.txt

# 按转录本坐标排序
sort -k1,1 -k2,2n condition_X_m5C_read_level.txt > condition_X_m5C_read_level_sorted.txt
sort -k1,1 -k2,2n condition_Y_m5C_read_level.txt > condition_Y_m5C_read_level_sorted.txt

# 模型 2 预测：汇总位点水平的 m5C 修饰概率与化学计量比（条件 X）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model2.py \
    -i condition_X_m5C_read_level_sorted.txt \
    -m ${MODEL_DIR}/CHEUI_m5C_model2.h5 \
    -c 0 \
    -d 0.3 0.7 \
    -n ${MIN_READS} \
    -o condition_X_m5C_site_level.txt

# 模型 2 预测：汇总位点水平的 m5C 修饰概率与化学计量比（条件 Y）
python3 ${CHEUI_PATH}/scripts/CHEUI_predict_model2.py \
    -i condition_Y_m5C_read_level_sorted.txt \
    -m ${MODEL_DIR}/CHEUI_m5C_model2.h5 \
    -c 0 \
    -d 0.3 0.7 \
    -n ${MIN_READS} \
    -o condition_Y_m5C_site_level.txt

# 生成差异分析配置文件，并使用 Mann-Whitney U 检验进行差异 m5C 分析
echo -e "X\tcondition_X_m5C_read_level_sorted.txt\nY\tcondition_Y_m5C_read_level_sorted.txt" \
    > m5C_diff_config.txt

python3 ${CHEUI_PATH}/scripts/CHEUI_diff.py \
    -i m5C_diff_config.txt \
    -n ${MIN_READS} \
    -o m5C_differential_results.txt
