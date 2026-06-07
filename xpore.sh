#!/bin/bash


# 参考转录组 FASTA 文件（cDNA，用于 minimap2 比对）
TRANSCRIPTOME="transcriptome.fa"

# 参考基因组注释文件（GTF 格式，推荐使用 GENCODE 或 ENSEMBL）
GTF="genome.gtf"

# 参考转录组 FASTA 文件（用于 xpore dataprep 坐标转换）
TRANSCRIPT_FASTA="transcriptome.fa"

# 并行线程数
THREADS=32

# 对照组样本目录（包含 fast5、fastq 文件）
CTRL_DIR="ctrl_rep1"

# 实验组样本目录
TREAT_DIR="treat_rep1"

# --------------------
# 【步骤 1】序列比对
# 使用 minimap2 将 fastq 比对到参考转录组
# -ax map-ont：Nanopore 长读长模式
# -uf：正链比对（直接 RNA 测序）
# --secondary=no：不输出次优比对结果
# --------------------

# 对照组比对
minimap2 \
    -ax map-ont \
    -uf \
    --secondary=no \
    -t ${THREADS} \
    ${TRANSCRIPTOME} \
    ${CTRL_DIR}/fastq/reads.fastq.gz \
    > ${CTRL_DIR}/aligned.sam

# 将 SAM 转换为排序后的 BAM 文件并建立索引
samtools view -Sb ${CTRL_DIR}/aligned.sam \
    | samtools sort -o ${CTRL_DIR}/aligned.sort.bam
samtools index ${CTRL_DIR}/aligned.sort.bam

# 实验组比对
minimap2 \
    -ax map-ont \
    -uf \
    --secondary=no \
    -t ${THREADS} \
    ${TRANSCRIPTOME} \
    ${TREAT_DIR}/fastq/reads.fastq.gz \
    > ${TREAT_DIR}/aligned.sam

samtools view -Sb ${TREAT_DIR}/aligned.sam \
    | samtools sort -o ${TREAT_DIR}/aligned.sort.bam
samtools index ${TREAT_DIR}/aligned.sort.bam

# --------------------
# 【步骤 2】nanopolish 索引与 eventalign
# 建立 fast5 与 fastq 的关联索引
# 将原始电信号比对到参考转录组，为 xpore 提供输入
# --signal-index：记录信号索引位置（xpore 必需）
# --scale-events：对电信号进行归一化缩放
# --------------------

# 对照组
nanopolish index \
    -d ${CTRL_DIR}/fast5 \              # fast5 原始信号目录
    ${CTRL_DIR}/fastq/reads.fastq.gz    # 对应的 fastq 文件

nanopolish eventalign \
    --reads ${CTRL_DIR}/fastq/reads.fastq.gz \
    --bam ${CTRL_DIR}/aligned.sort.bam \
    --genome ${TRANSCRIPTOME} \         # 参考转录组文件
    --signal-index \
    --scale-events \
    --summary ${CTRL_DIR}/summary.txt \
    --threads ${THREADS} \
    > ${CTRL_DIR}/eventalign.txt

# 实验组
nanopolish index \
    -d ${TREAT_DIR}/fast5 \
    ${TREAT_DIR}/fastq/reads.fastq.gz

nanopolish eventalign \
    --reads ${TREAT_DIR}/fastq/reads.fastq.gz \
    --bam ${TREAT_DIR}/aligned.sort.bam \
    --genome ${TRANSCRIPTOME} \
    --signal-index \
    --scale-events \
    --summary ${TREAT_DIR}/summary.txt \
    --threads ${THREADS} \
    > ${TREAT_DIR}/eventalign.txt

# --------------------
# 【步骤 3】xpore dataprep 数据预处理
# 将 nanopolish eventalign 输出转换为 xpore 所需的 JSON 格式
# --gtf_or_gff：基因组注释文件，用于转录组坐标到基因组坐标的映射
# --transcript_fasta：转录本序列文件
# --genome：启用基因组坐标模式输出
# --n_processes：并行进程数
# --------------------

# 对照组预处理
xpore dataprep \
    --eventalign ${CTRL_DIR}/eventalign.txt \   # 输入 eventalign 文件
    --gtf_or_gff ${GTF} \                       # 基因组注释 GTF 文件
    --transcript_fasta ${TRANSCRIPT_FASTA} \    # 转录本 FASTA 文件
    --out_dir ${CTRL_DIR}/dataprep \            # 输出目录
    --genome \                                  # 输出基因组坐标
    --n_processes ${THREADS}

# 实验组预处理
xpore dataprep \
    --eventalign ${TREAT_DIR}/eventalign.txt \
    --gtf_or_gff ${GTF} \
    --transcript_fasta ${TRANSCRIPT_FASTA} \
    --out_dir ${TREAT_DIR}/dataprep \
    --genome \
    --n_processes ${THREADS}

# --------------------
# 【步骤 4】编写 YAML 配置文件
# 指定实验设计、各样本数据目录、输出目录及统计方法参数
# --------------------

cat > config.yml << EOF
data:
    Control:                                    # 对照组名称
        rep1: ./${CTRL_DIR}/dataprep            # 对照组 dataprep 输出目录
    Treatment:                                  # 实验组名称
        rep1: ./${TREAT_DIR}/dataprep           # 实验组 dataprep 输出目录

out: ./xpore_output                             # 差异修饰分析结果输出目录

method:
    prefiltering:                               # 预过滤以排除不太可能存在差异修饰的位点
        method: t-test                          # 使用 t 检验进行预过滤
        threshold: 0.1                          # 预过滤 p-value 阈值
EOF

# --------------------
# 【步骤 5】xpore diffmod 差异修饰分析
# 基于贝叶斯高斯混合模型对每个位点进行差异修饰检测
# 输出 diffmod.table：包含所有检测位点的差异修饰率和 p-value
# --n_processes：并行进程数
# --save_models：保存每个位点的模型参数（可选）
# --------------------

mkdir -p xpore_output

xpore diffmod \
    --config config.yml \           # 输入 YAML 配置文件
    --n_processes ${THREADS} \      # 并行进程数
    --save_models                   # 保存每个位点的模型参数

# --------------------
# 【步骤 6】结果后处理（可选）
# 对 diffmod.table 结果进行整理，添加多重检验校正等信息
# --------------------

xpore postprocessing \
    --diffmod_dir xpore_output      # 输入 diffmod 结果目录
