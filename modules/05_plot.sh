#!/bin/bash
# modules/05_plot.sh

BAM_FILE="${DIR_BAM}/${FOLDER_BASENAME}_aligned.sorted.bam"
PY_SCRIPT="${PY_SCRIPT_DIR}/WGSmapping.py"

if [[ ! -f "$BAM_FILE" ]]; then echo "[ERROR] BAM missing"; exit 1; fi

echo "[MODULE: PLOT] Running Linear Pile-up visualization..."
echo "  - Window Radius: +/- $WINDOW_SIZE bp"

# === [修改点] 增加 --window 参数传递 ===
python "$PY_SCRIPT" \
    --bam "$BAM_FILE" \
    --chromosome "$CHROMOSOME" \
    --center "$CENTER_POSITION" \
    --step "$STEP_SIZE" \
    --window "$WINDOW_SIZE" \
    --background "$BACKGROUND_ANALYSIS" \
    --output "$DIR_RESULT" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then exit 1; fi