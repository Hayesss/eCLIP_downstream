# eCLIP_downstream

基于 **hg38 / GRCh38 + GENCODE v38** 的 SerRS/SARS1 eCLIP 结合区间序列提取与注释版本核查。

## 功能

- 下载 GENCODE v38 官方参考及基因/转录本 ID 元数据，核对官方 MD5。
- 按原始位点提取参考序列，生成总长 500 bp 和 1000 bp 的扩展窗口。
- 根据正负链输出 5′→3′ 方向序列；负链反向互补，使用 DNA 字母。
- 分别输出 FASTA、TSV，保留基因/转录本、坐标、结合区间在序列中的偏移和质控信息。
- 比较带版本 ID 与多个 GENCODE release 的匹配率。

## 安装

需要 Python 3.10 或更新版本，建议 64 位 Python。已在 Python 3.12、openpyxl 3.1.5 验证。

```bash
python -m pip install -r requirements.txt
```

## 使用

### 1. 准备参考文件

```bash
python scripts/download_reference.py --workers 8
```

默认缓存于 `references/gencode_v38/`。参考压缩文件约 845 MB，解压后约 3.15 GB，准备至少 5 GB 可用空间。分段下载可以重试并复用已完成分段；完成后核对官方 MD5。元数据会一并下载，不需要事先准备目录或校验文件。

已有缓存时可指定位置：

```bash
python scripts/download_reference.py --output-dir /path/to/reference-cache
```

只需要 ID 元数据时可加 `--metadata-only`。

### 2. 提取三套序列

```bash
python scripts/extract_sequences.py --input "data/sars1_eclip_binding_sit4es.xlsx" --expected-count 50613
```

默认工作表为 `Sheet1`。官方补充工作簿可指定 `--sheet "Table S3"`。表头必须在第 2 行，数据从第 3 行开始，采用 Table S3 的列结构。

可以指定参考目录及输出目录：

```bash
python scripts/extract_sequences.py --input "data/sites.xlsx" --reference-dir "/path/to/reference-cache" --output-dir "results/my_sequences"
```

`--expected-count` 是可选的完整性检查。文件保留原始顺序，不合并重复位点、不增加显著性筛选。必需列包括 `chr/start/end/name/strand/gc`、基因和转录本注释、区域注释、`q_min/q_max`、`enrichment_l2or_mean`、`input_sum/clip_sum`。ID 必须与 v38 官方清单兼容；原表 GC 一致率低于 99% 会停止导出。

默认输出位于 `results/sequences_hg38_gencodev38/`：

| 输出前缀 | 内容 |
| --- | --- |
| `SARS1_eCLIP_binding_sites_hg38_GENCODEv38` | 原始结合区间，FASTA 和 TSV |
| `SARS1_eCLIP_centered_500bp_hg38_GENCODEv38` | 总长 500 bp，FASTA 和 TSV |
| `SARS1_eCLIP_centered_1000bp_hg38_GENCODEv38` | 总长 1000 bp，FASTA 和 TSV |
| `sequence_extraction_QC.json` | 参考来源、GC/长度检查、ID 匹配及 SHA-256 |
| `SARS1_eCLIP_sequences_hg38_GENCODEv38_all.zip` | 六个数据文件及说明、质控文件 |

重复运行会覆盖所选输出目录中的同名生成文件；不同分析应指定不同 `--output-dir`。请勿使用 Python `-O`，以保留科学验证断言。

### 3. 可选：重现注释版本比对

```bash
python scripts/audit_gencode_versions.py --input "data/sars1_eclip_binding_sit4es.xlsx"
```

默认比较 v29 和 v37–v43，可用 `--releases 37 38 39 40` 调整。公开 ID 列表下载后在本地与输入比较，输入文件、ID 和坐标不发送至外部服务。输出 TSV 和 JSON 位于 `results/reference_audit/`。

## 序列与坐标约定

- 使用 **0-based、左闭右开** 的基因组坐标，区间长度为 `end-start`。
- 500/1000 bp 是包含结合区间在内的最终总长度。两侧尽量均分；奇数时额外 1 bp 位于 RNA 方向的 3′侧。
- 扩展完整保留原始结合区间。线性染色体边缘越界时向内平移，chrM 越界时按环状序列拼接，均有标记。
- 序列为连续基因组区间，可能含内含子；未按成熟转录本拼接、未加入样本特异变异。
- 输出的是 eCLIP 结合区间序列，不表示最小结合基序或每个碱基均直接接触 SerRS。
- 使用 `sequence[peak_offset0_in_sequence:peak_end_offset0_in_sequence]` 可取回原始结合序列。
- 少量参考未知碱基以 N 保留，并记录 `N_count`。
- 现有结合位点注释被保留并核对 v38 ID；不重新计算外显子/UTR 注释，不依赖原表 uga_* 列。
- TSV 可通过 Excel“数据→从文本/CSV”导入，分隔符选制表符。基因名和 ID 列按文本读取，避免自动日期转换。

## 验证

```bash
python -m unittest discover -s tests -v
```

测试使用合成 DNA，不需要网络或用户数据；覆盖 FASTA 索引、正负链、窗口长度、染色体边缘、环状 chrM 和分段下载响应验证。全量提取还会逐条核对原始序列嵌入关系以及 FASTA/TSV 一致性。

[参考版本依据](docs/reference_versions.md)记录原研究来源和版本核查结果。原始表、参考基因组、下载缓存及生成的大文件不纳入 Git。
