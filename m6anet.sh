#!/bin/bash
# ============================================================
# m6A 修饰位点检测流程（nanopolish + m6anet）
# ============================================================

# --------------------
# 【步骤 1】nanopolish 索引
# 将原始 fast5 信号文件与 fastq 文件关联，建立索引
# --------------------
echo ">>> 步骤 1：建立 nanopolish 索引..."

nanopolish index \
    -d /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-2/20230818-NPL231153-P4-PAQ94845/PAQ94845/fast5_pass \
    # fast5 原始信号文件目录
    /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-2/20230818-NPL231153-P4-PAQ94845/PAQ94845/T-2.fastq
    # 对应的 fastq 文件

echo "    索引建立完成！"

# --------------------
# 【步骤 2】nanopolish eventalign
# 将电信号与参考转录组进行比对，输出每个碱基对应的电信号特征
# --------------------
echo ""
echo ">>> 步骤 2：nanopolish eventalign 电信号比对..."

nanopolish eventalign \
    --reads /data1/zhouyiming/qian_3/WHXWZB-2023080124A/raw_data/Nanopore/T-1/20230818-NPL231152-P4-PAQ94711/PAQ94711/T-1.fastq \
    # 输入 fastq 文件
    --bam /data1/zhouyiming/qian_3/drsseqres/transcript/fq/T-1_fq/T-1.bam \
    # 比对到参考转录组的 BAM 文件
    --genome /data2/backup/share_1423_backup/liuqi/ribo/ribotoolkit/db/mRNA/osa_IRGSP_1.txdb.fa \
    # 参考转录组 FASTA 文件
    --scale-events \
    # 对电信号进行归一化缩放
    --summary T-1_summary.txt \
    # 输出每条 read 的处理摘要文件
    --signal-index \
    # 在输出中记录信号索引（m6anet 必需）
    --threads 50 \
    # 使用 50 个线程并行运行
    > T-1_eventalign.txt
    # 输出 eventalign 结果文件

echo "    eventalign 比对完成！结果保存至 T-1_eventalign.txt"

# --------------------
# 【步骤 3】m6anet 数据预处理
# 将 eventalign 输出文件转换为 m6anet 推断所需的格式
# --------------------
echo ""
echo ">>> 步骤 3：m6anet dataprep 数据预处理..."

m6anet dataprep \
    --eventalign /data1/zhouyiming/DT4/DT4-2-1/DT4-2_eventalign.txt \
    # 输入 nanopolish eventalign 结果文件
    --out_dir dataprep \
    # 预处理结果输出目录
    --n_processes 8 \
    # 使用 8 个进程并行处理
    --readcount_max 2000000
    # 每个位点最多读取的 read 数量上限

echo "    数据预处理完成！结果保存至 dataprep/ 目录"

# --------------------
# 【步骤 4】m6anet 推断
# 基于预训练模型对每个位点进行 m6A 修饰概率预测
# --------------------
echo ""
echo ">>> 步骤 4：m6anet inference m6A 位点预测..."

m6anet inference \
    --input_dir dataprep \
    # 输入目录（步骤 3 的输出目录）
    --out_dir run \
    # 预测结果输出目录
    --pretrained_model arabidopsis_RNA002 \
    # 使用拟南芥 RNA002 预训练模型
    --n_processes 8 \
    # 使用 8 个进程并行推断
    --num_iterations 1000
    # 模型推断迭代次数

echo "    m6A 位点预测完成！结果保存至 run/ 目录"

# --------------------
# 【完成】
# --------------------
echo ""
echo "============================================"
echo "  全流程运行完毕！"
echo "  最终预测结果：run/"
echo "============================================"
