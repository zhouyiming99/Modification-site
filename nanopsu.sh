#!/bin/bash


# 原始 fast5 数据目录（未碱基识别时使用）
RAW_FAST5="fast5"

# 碱基识别输出目录
FASTQ_DIR="fastq"

# 参考基因组文件
REFERENCE="reference.fa"

# Nanopore 使用的流动池型号
FLOWCELL="FLO-MIN106"

# Nanopore 直接 RNA 测序试剂盒版本
KIT="SQK-RNA002"

# 每个碱基识别进程的 CPU 线程数
CPU_THREADS=4

# 碱基识别并行进程数
NUM_CALLERS=4


echo ">>> 检查 nanopsu 是否可用..."
if ! command -v nanopsu &> /dev/null; then
    echo "[错误] 未找到 nanopsu 命令，请先安装该工具："
    echo "       cd <Nanopore_psU目录> && pip install ."
    exit 1
fi
echo "    nanopsu 已安装，继续执行..."

# 【步骤 1】碱基识别（Basecalling）
# 如果数据在测序时已进行碱基识别，可跳过此步骤
# 使用 guppy_basecaller 将原始 fast5 信号转换为 fastq 序列

echo ""
echo ">>> 步骤 1：碱基识别（Basecalling）"
echo "    注意：若数据已完成碱基识别，请注释掉此步骤"

guppy_basecaller \
    --input_path ${RAW_FAST5} \      # 输入原始 fast5 文件目录
    --recursive \                     # 递归搜索子目录中的 fast5 文件
    --save_path ${FASTQ_DIR} \        # 输出 fastq 文件的保存目录
    --records_per_fastq 0 \           # 0 表示所有 reads 写入同一个 fastq 文件
    --flowcell ${FLOWCELL} \          # 指定流动池型号
    --kit ${KIT} \                    # 指定测序试剂盒
    --qscore_filtering \              # 开启质量分数过滤
    --min_qscore 7 \                  # 最低质量分数阈值设为 7
    --cpu_threads_per_caller ${CPU_THREADS} \   # 每个识别进程的 CPU 线程数
    --num_callers ${NUM_CALLERS}      # 并行运行的识别进程数

echo "    碱基识别完成，fastq 文件已保存至：${FASTQ_DIR}"


# 【步骤 2】序列比对与 pileup（Alignment）
# 将 fastq reads 比对到参考基因组，并堆叠 reads

echo ""
echo ">>> 步骤 2：序列比对与 pileup（Alignment）"
echo "    

nanopsu alignment \
    -i ${FASTQ_DIR} \       # 输入 fastq 文件目录（fastq 文件需直接位于该目录下）
    -r ${REFERENCE}          # 参考基因组 FASTA 文件

echo "    比对完成，结果保存至 alignment/ 目录"
echo "    包含子目录：alignment/plus_strand/ 和 alignment/minus_strand/"


# 【步骤 3】去除内含子跳跃区域（Remove Intron Gaps）
# 由于 samtools 设计，mpileup 文件中剪接 reads 跳跃区域会被填充 '>' 或 '<'
# 这些并非真实碱基，需要去除以避免影响后续分析

echo ""
echo ">>> 步骤 3：去除 mpileup 文件中的内含子跳跃区域"

nanopsu rm_intron \
    -i alignment/           # 输入 alignment 目录路径

echo "    完成！在 plus_strand/ 和 minus_strand/ 目录下生成 collect_pile_no_intron.txt"

# 【步骤 4】提取特征（Feature Extraction）
# 提取所有 U 位点的特征信息
# 覆盖度阈值为 20 reads，只有覆盖深度 >20 的 U 位点才会被分析

echo ""
echo ">>> 步骤 4：提取 U 位点特征（Feature Extraction）"


nanopsu get_features \
    -i alignment/           # 输入 alignment 目录路径

echo "    特征提取完成！结果保存至 alignment/features.csv"
echo "    该文件包含正链和负链上所有符合条件的 U 位点信息"

# --------------------
# 【步骤 5】预测假尿嘧啶（psU Prediction）
# 基于机器学习模型预测每个 U 位点的假尿嘧啶（Ψ）概率
# 输出文件 prediction.csv 包含：参考链、位置、碱基类型、覆盖度、U概率、psU概率
# --------------------
echo ""
echo ">>> 步骤 5：假尿嘧啶位点预测（psU Prediction）"

nanopsu predict \
    -i alignment/features.csv   # 输入特征文件（由步骤 4 生成）

echo "    预测完成！结果保存至 alignment/prediction.csv"
echo ""
echo "    prediction.csv 每行包含以下字段："
echo "    参考链 | 位置 | 碱基类型 | 覆盖度 | U概率 | psU概率"

# --------------------
# 【完成】
# --------------------
echo ""
echo "============================================"
echo "  NanoPsu 全流程运行完毕！"
echo "  最终预测结果：alignment/prediction.csv"
echo "============================================"
