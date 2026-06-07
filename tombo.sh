#!/bin/bash

# 依赖环境：
# h5py < 3
# numpy < 1.20
# scipy
# cython
# setuptools >= 18.0
# mappy >= 2.10
# future
# tqdm

# 参考转录组路径
ref=path_to_reference_transcriptome

# basecalled fast5文件路径
fast5=path_to_basecalled_fast5_files

# 重新比对原始信号（resquiggling）
tombo resquiggle $fast5/workspace $ref \
    --rna \
    --corrected-group RawGenomeCorrected_000 \
    --basecall-group Basecall_1D_000 \
    --overwrite \
    --processes 16 \
    --fit-global-scale \
    --include-event-stdev

# 从头检测RNA修饰位点
tombo detect_modifications de_novo \
    --fast5-basedirs $fast5/workspace \
    --statistics-file-basename sample \
    --corrected-group RawGenomeCorrected_000 \
    --processes 16

# 输出RRACH基序下m6A位点的统计结果（覆盖度、校正分数、修饰比例）
tombo text_output browser_files \
    --fast5-basedirs $fast5/workspace \
    --statistics-filename sample.tombo.stats \
    --browser-file-basename sample_rrach \
    --genome-fasta $ref \
    --motif-descriptions RRACH:3:m6A \
    --file-types coverage dampened_fraction fraction \
    --corrected-group RawGenomeCorrected_000
