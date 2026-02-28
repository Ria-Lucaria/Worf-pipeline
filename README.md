# WORF-Seq Pipeline

WORF-Seq 是一个基于 shell + Python 的全基因组测序（WGS）小型流水线，用于：
- 质控：`fastp`
- 比对：`minimap2`
- BAM 处理：`samtools sort/index`
- 可视化与统计：`scripts/WGSmapping.py`（Matplotlib + NumPy）

该流水线按模块顺序运行：`01_qc.sh → 02_align.sh → 04_process.sh → 03_db_ingest.sh → 05_plot.sh`，由入口脚本 `worf.sh` 调度。

## 目录结构
```
worf.sh
modules/
  01_qc.sh
  02_align.sh
  03_db_ingest.sh
  04_process.sh
  05_plot.sh
references/
scripts/
  WGSmapping.py
  bam2mongo.py
  query_worf.py
  fetch_genes.py
  bam2mongo.py
server/
  main.py
  database.py
  models.py
  README.md
worf_tools.sh
```

## 环境准备
建议使用 Conda 创建并激活环境（默认环境名：`worf_env`）。

```bash
conda env create -f environment.yml
conda activate worf_env
```

如果只需要 Python 包，也可以使用 pip：

```bash
python -m pip install -r requirements.txt
```

## 外部工具依赖
- `fastp`（质控）
- `minimap2`（比对）
- `samtools`（BAM 处理与区域计数）
- `pysam`（Python 读写 BAM，入库脚本使用）

这些工具已在 [environment.yml](environment.yml) 中列出。若使用 Conda 手动安装：

```bash
conda install -c bioconda fastp minimap2 samtools pysam
```

注意：运行时需要确保相关工具在 `PATH` 中或已激活包含这些工具的环境。

## 运行方式（Pipeline）
入口脚本：`worf.sh`

```bash
bash worf.sh -f <input_folder> -c <chrom> -p <center_bp> \
  [-w window_bp] [-s step_bp] [-b true|false] [-o outdir] [-r]
```

批处理模式：

```bash
bash worf.sh -B batch.csv [-w window_bp] [-s step_bp] [-b true|false] [-o batch_outdir]
```

- `-f`：输入目录（包含原始 FASTQ）
- `-c`：染色体（例如 `chr6` 或 `6`，需与 BAM 头的 `@SQ SN:` 一致）
- `-p`：中心位置（bp）
- `-w`：目标窗口半径（bp，默认 `10000`）
- `-s`：步长（bp，默认 `100000`，用于全染色体背景统计）
- `-b`：是否进行背景分析（默认 `true`）
- `-o`：输出目录（默认 `<input_folder>/output`）
- `-r`：强制重跑（忽略检查点）

示例：

```bash
bash worf.sh -f /data/.../UDI001 -c chr6 -p 31429000 -w 10000
```

批处理示例：

```bash
cp batch.csv.example batch.csv
bash worf.sh -B batch.csv
```

`batch.csv` 列顺序：

```text
task_id,input_folder,chromosome,center,window,step,background,output_dir,force_run
```

说明：
- `task_id`、`window`、`step`、`background`、`output_dir`、`force_run` 可留空（使用默认值）。
- 若 `output_dir` 为空，默认输出到 `<input_folder>/output/<task_id>`。
- 批处理结束后会在项目根目录生成 `batch_run_summary_*.csv`。

### Checkpoints 与重复绘图
- 前三步（`01_qc.sh`, `02_align.sh`, `04_process.sh`）使用 checkpoint 机制（`.01_qc.done` 等）。
- `05_plot.sh` 不再创建/检查 checkpoint，允许使用不同 `-c/-p/-w` 参数反复绘图。

### 绘图输出命名
- 由 `scripts/WGSmapping.py` 生成的图像文件名包含 `-c`（染色体）、`-p`（中心位置）、`-w`（窗口半径）与时间戳，便于区分多次绘图。

### 每任务分析输出（新增）
- 命中 reads 明细：`*_hit_reads.csv`（每条命中 read 的位置、长度、MAPQ、CIGAR、链方向、overlap 信息）。
- 任务 summary：`*_summary.txt`（包含 `total_hit_reads`、`unique_hit_reads`、`mean_depth`、`max_depth`、覆盖比例等统计）。

## 数据入库与查询
- 第三步 `modules/03_db_ingest.sh` 调用 `scripts/bam2mongo.py` 将 BAM 映射元数据写入 MongoDB（默认 `127.0.0.1:30001`）。
- 查询与可视化结果可通过服务端 API 访问，详见 [server/README.md](server/README.md)。

### 独立工具脚本：`worf_tools.sh`
- 功能：提供交互式/命令行方式进行数据库入库与查询。
- 依赖：`scripts/bam2mongo.py`、`scripts/query_worf.py`，以及 MongoDB 服务（默认 `127.0.0.1:30001`）。
- 用法：
  - 打开交互式菜单：`./worf_tools.sh`
  - 手动入库：`./worf_tools.sh ingest <bam_path> [sample_id]`
  - 交互式查询（TUI）：`./worf_tools.sh query`

## 服务端 API（FastAPI）
启动服务：

```bash
uvicorn server.main:app --reload --host 0.0.0.0 --port 8018
```

文档与测试：
- Swagger UI: `http://localhost:8018/docs`
- 端点：健康检查 `/`、样品列表 `/samples`、区域查询 `/query/{sample_id}`（参数：`chrom`、`start`、`end`、`limit`）。

## 可视化脚本说明（`scripts/WGSmapping.py`）
- 使用 `samtools` 统计区域内 Reads 并生成覆盖度图与堆叠图。
- 性能优化：优先 `BrokenBarHCollection`，兼容回退 `PatchCollection` + `Rectangle`。

## 常见问题
- 找不到 `fastp/minimap2/samtools/pysam`：请确保 Conda 环境已激活，或手动安装这些工具。
- 找不到原始 FASTQ：确认文件名是否满足约定（`<basename>_raw_1.fq.gz` 等）。
- 参考缺失：在 `modules/02_align.sh` 中修改 `REF_DIR`，指向你的参考文件位置。
- 图像为空或覆盖度为零：确认 BAM 已正确生成并索引，染色体名称与 BAM 头一致（`@SQ SN:<chrom>`）。

## 致谢
- fastp, minimap2, samtools, pysam; Matplotlib 与 NumPy 用于绘图与数值计算。