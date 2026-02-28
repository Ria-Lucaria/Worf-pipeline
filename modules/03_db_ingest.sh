#!/bin/bash
set -e  # 遇到错误立即退出
set -o pipefail

# ===============================================================
# Module 03: Database Ingestion
# Description: 将比对后的 BAM 文件元数据导入 MongoDB
# Author: Guo Zehua (Worf Pipeline)
# ===============================================================

# 接收参数 (由 master script 传入)
SAMPLE_ID=$1
CONFIG_FILE=$2  # 预留，如果未来有配置文件
DATE_TAG=$(date +%Y%m%d) # 自动生成今天的日期标签

# ================= 配置区域 =================
# 输入目录 (必须与 04_process.sh 的输出保持一致)
INPUT_BAM_DIR="${DIR_BAM:-results/03_bam}"

# Python 脚本位置
INGEST_SCRIPT="${PY_SCRIPT_DIR}/bam2mongo.py"

# MongoDB 配置
MONGO_PORT=30001
MONGO_HOST="127.0.0.1"

# 构造文件名 (根据你的命名习惯)
BAM_FILE="${INPUT_BAM_DIR}/${SAMPLE_ID}_aligned.sorted.bam"

# 构造数据库中的 Sample Name (建议: 日期-样品名，方便区分批次)
DB_SAMPLE_NAME="${DATE_TAG}-${SAMPLE_ID}"
# ===========================================

# 1. 检查输入文件
if [ ! -f "$BAM_FILE" ]; then
    echo "[Error] BAM file not found: $BAM_FILE"
    echo "Please check if Module 02 finished successfully."
    exit 1
fi

echo "=================================================="
echo "[Step 03] Starting Database Ingestion"
echo "Target Sample: ${SAMPLE_ID}"
echo "DB Entry Name: ${DB_SAMPLE_NAME}"
echo "BAM Source   : ${BAM_FILE}"
echo "Start Time   : $(date)"
echo "=================================================="

# 2. 调用 Python 入库脚本
# 注意：这里假设你在 worf.sh 中已经激活了 conda 环境
python "$INGEST_SCRIPT" \
    --bam "$BAM_FILE" \
    --sample "$DB_SAMPLE_NAME" \
    --port "$MONGO_PORT" \
    --host "$MONGO_HOST"

# 3. 结果检查
if [ $? -eq 0 ]; then
    echo "--------------------------------------------------"
    echo "[Success] Ingestion complete for ${SAMPLE_ID}"
    echo "End Time     : $(date)"
    echo "=================================================="
else
    echo "[Error] Python script failed."
    exit 1
fi