# 参考版本依据

本项目固定使用 **hg38 / GRCh38 坐标和 GENCODE v38 注释**。

原研究：Liu et al., *Nucleic Acids Research* (2023), DOI [10.1093/nar/gkad773](https://doi.org/10.1093/nar/gkad773)。补充 Table S3 包含 50,613 个 SerRS eCLIP 结合位点。原文提供 hg38 的 Skipper 分析轨道，但 eCLIP 方法段未明确写出 GENCODE release 或 GRCh38 补丁号。正文的 human transcriptome (v19) 属于 ribosome profiling，不能套用于 eCLIP。

我们对表中完整的带版本 ID 进行了本地精确匹配。以下结果为版本兼容性证据，并非作者对原始运行配置的直接确认。

| GENCODE release | gene_id（7,166）匹配 | transcript_ids（24,780）匹配 | uga_gene_id（5,674）匹配 |
| --- | ---: | ---: | ---: |
| 29 | 820 | 17,713 | 173 |
| 37 | 6,839 | 24,436 | 4,743 |
| 38 | 7,166 | 24,780 | 5,012 |
| 39 | 6,375 | 24,285 | 5,581 |
| 40 | 6,283 | 24,178 | 5,674 |
| 41 | 6,034 | 23,960 | 5,440 |
| 42 | 5,869 | 23,884 | 5,340 |
| 43 | 5,817 | 23,843 | 5,287 |

结合位点注释高度支持 v38，UGA 注释高度支持 v40 或对应的派生注释。本项目根据结合位点坐标提取序列，不使用 uga_* 列决定序列位置。ID 一致不等于所有外显子、UTR 或 stop_codon 坐标已经逐一核查。

GENCODE v38 官方发布标签为 GRCh38.p13；本项目明确选用其中的 **GRCh38 primary assembly** FASTA，而非包含全部 patch/alternate loci 的 FASTA。这不能证明作者原始比对 FASTA 的补丁版本。

固定参考文件：`GRCh38.primary_assembly.genome.fa.gz`

- 官方发布：https://www.gencodegenes.org/human/release_38.html
- 官方下载：https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/
- 压缩文件 MD5：`71b78189dc2a7b7e7af802949f168aa7`
- 压缩文件 SHA-256：`0c792010345db1834a29e2df9bf053ea30d3b5ebab3536710632374ef73f7f6a`
- GRC 关于补丁不改变主染色体坐标的说明：https://www.ncbi.nlm.nih.gov/grc/help/patches/

完整输入表的实测验证：50,613/50,613 个结合区间的 GC 含量在原表精度内一致（容差 1e-6）；0-based、右端不包含的坐标解释得到支持。原始序列长 1–100 bp，扩展结果严格为 500 bp 和 1000 bp。所有序列按源表 strand 输出，负链反向互补。上述数值是指定数据的核查结果，并非对任意输入的保证。
