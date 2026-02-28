#!/bin/bash
# worf_pipeline/worf.sh
umask 000
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${SCRIPT_DIR}/modules"
PY_SCRIPT_DIR="${SCRIPT_DIR}/scripts"

# === 默认参数 ===
STEP_SIZE=100000
WINDOW_SIZE=10000  # 默认精细窗口半径 (10kb)
BACKGROUND_ANALYSIS=true
FORCE_RUN=false
BATCH_CSV=""

usage() {
    echo "Usage:"
    echo "  Single task: $0 -f <folder> -c <chrom> -p <center> [-w window] [-s step] [-b background] [-o outdir] [-r]"
    echo "  Batch mode : $0 -B <batch.csv> [-w window] [-s step] [-b background] [-o outdir]"
    echo ""
    echo "Batch CSV columns:"
    echo "  task_id,input_folder,chromosome,center,window,step,background,output_dir,force_run"
    echo "  -w  Target window radius (bp) for linear pile-up plot (default: 10000)"
    echo "  -r  Force re-run (overwrite existing checkpoints)"
    exit 1
}

# === 参数解析 ===
while getopts ":f:c:p:w:s:b:o:rB:" opt; do
    case $opt in
        f) INPUT_FOLDER="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        p) CENTER_POSITION="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        s) STEP_SIZE="$OPTARG" ;;
        b) BACKGROUND_ANALYSIS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) FORCE_RUN=true ;;
        B) BATCH_CSV="$OPTARG" ;;
        \?) usage ;;
        :) usage ;;
    esac
done

GLOBAL_OUTPUT_DIR="$OUTPUT_DIR"

trim() {
    local value="$1"
    value="${value#${value%%[![:space:]]*}}"
    value="${value%${value##*[![:space:]]}}"
    echo "$value"
}

is_truthy() {
    local value
    value=$(echo "$1" | tr '[:upper:]' '[:lower:]')
    [[ "$value" == "1" || "$value" == "true" || "$value" == "yes" || "$value" == "y" ]]
}

setup_task_context() {
    INPUT_FOLDER="$1"
    CHROMOSOME="$2"
    CENTER_POSITION="$3"
    WINDOW_SIZE="$4"
    STEP_SIZE="$5"
    BACKGROUND_ANALYSIS="$6"
    OUTPUT_DIR="$7"
    FORCE_RUN="$8"

    FOLDER_BASENAME=$(basename "$INPUT_FOLDER")
    SAMPLE_ID="${FOLDER_BASENAME}"

    export DIR_QC="${OUTPUT_DIR}/01_qc"
    export DIR_ALIGN="${OUTPUT_DIR}/02_align"
    export DIR_BAM="${OUTPUT_DIR}/03_bam"
    export DIR_RESULT="${OUTPUT_DIR}/04_result"

    if ! mkdir -p "$OUTPUT_DIR" "$DIR_QC" "$DIR_ALIGN" "$DIR_BAM" "$DIR_RESULT" 2>/dev/null; then
        echo "[ERROR] Cannot create output directories under: $OUTPUT_DIR" >&2
        return 1
    fi

    LOG_FILE="${OUTPUT_DIR}/pipeline.log"

    export INPUT_FOLDER OUTPUT_DIR LOG_FILE SCRIPT_DIR PY_SCRIPT_DIR
    export CHROMOSOME CENTER_POSITION STEP_SIZE WINDOW_SIZE BACKGROUND_ANALYSIS
    export FOLDER_BASENAME SAMPLE_ID FORCE_RUN
}

log_info() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
    if [[ -n "$LOG_FILE" && -d "$(dirname "$LOG_FILE")" ]]; then
        echo "$msg" | tee -a "$LOG_FILE"
    else
        echo "$msg"
    fi
}

log_error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1"
    if [[ -n "$LOG_FILE" && -d "$(dirname "$LOG_FILE")" ]]; then
        echo "$msg" | tee -a "$LOG_FILE"
    else
        echo "$msg"
    fi
}

# === 模块执行函数 ===
run_module() {
    local step_name="$1"
    local script_path="${MODULE_DIR}/${step_name}"
    local checkpoint_file="${OUTPUT_DIR}/.${step_name%.*}.done"

    # 检查脚本是否存在
    if [[ ! -f "$script_path" ]]; then
        log_error "Module script not found: $script_path"
        return 1
    fi

    # Checkpoint 检查 (05_plot 总是运行，不跳过)
    if [[ "$step_name" != "05_plot.sh" && "$FORCE_RUN" == "false" && -f "$checkpoint_file" ]]; then
        log_info "Skipping $step_name (Checkpointed)"
        return 0
    fi

    log_info "Running $step_name for Sample: $SAMPLE_ID..."
    
    # [修改点] 执行时传入 SAMPLE_ID 作为 $1 参数
    # 这样 modules/03_db_ingest.sh 里的 SAMPLE_ID=$1 就能正确获取
    if /bin/bash "$script_path" "$SAMPLE_ID"; then
        log_info "$step_name Completed."
        # 成功后创建 checkpoint
        if [[ "$step_name" != "05_plot.sh" ]]; then
            touch "$checkpoint_file"
        fi
    else
        log_error "$step_name Failed."
        return 1
    fi
}

run_single_task() {
    local input_folder="$1"
    local chromosome="$2"
    local center_position="$3"
    local window_size="$4"
    local step_size="$5"
    local background_analysis="$6"
    local output_dir="$7"
    local force_run="$8"

    if ! setup_task_context "$input_folder" "$chromosome" "$center_position" "$window_size" "$step_size" "$background_analysis" "$output_dir" "$force_run"; then
        return 1
    fi

    log_info "Task Start: sample=${SAMPLE_ID}, region=${CHROMOSOME}:${CENTER_POSITION}, output=${OUTPUT_DIR}"

    run_module "01_qc.sh" || return 1
    run_module "02_align.sh" || return 1
    run_module "04_process.sh" || return 1
    run_module "03_db_ingest.sh" || return 1
    run_module "05_plot.sh" || return 1

    chmod -R 777 "$OUTPUT_DIR" 2>/dev/null || true
    log_info "Task Done."
}

run_batch() {
    local batch_file="$1"
    local batch_ts
    local batch_summary
    batch_ts=$(date +%Y%m%d_%H%M%S)
    batch_summary="${SCRIPT_DIR}/batch_run_summary_${batch_ts}.csv"

    if [[ ! -f "$batch_file" ]]; then
        echo "[ERROR] Batch CSV not found: $batch_file"
        exit 1
    fi

    echo "task_id,input_folder,chromosome,center,window,step,background,output_dir,force_run,status,message" > "$batch_summary"

    local line_no=0
    local success_count=0
    local fail_count=0

    while IFS=',' read -r task_id input_folder chromosome center window step background output_dir force_run; do
        line_no=$((line_no + 1))

        # 跳过空行和注释
        [[ -z "$task_id$input_folder$chromosome$center$window$step$background$output_dir$force_run" ]] && continue
        [[ "$(trim "$task_id")" =~ ^# ]] && continue

        # 跳过表头
        local normalized_header
        normalized_header=$(echo "$(trim "$task_id"),$(trim "$input_folder"),$(trim "$chromosome"),$(trim "$center")" | tr '[:upper:]' '[:lower:]')
        if [[ "$normalized_header" == "task_id,input_folder,chromosome,center" || "$normalized_header" == "input_folder,chromosome,center," ]]; then
            continue
        fi

        task_id=$(trim "$task_id")
        input_folder=$(trim "$input_folder")
        chromosome=$(trim "$chromosome")
        center=$(trim "$center")
        window=$(trim "$window")
        step=$(trim "$step")
        background=$(trim "$background")
        output_dir=$(trim "$output_dir")
        force_run=$(trim "$force_run")

        [[ -z "$task_id" ]] && task_id="task_${line_no}"
        [[ -z "$window" ]] && window="$WINDOW_SIZE"
        [[ -z "$step" ]] && step="$STEP_SIZE"
        [[ -z "$background" ]] && background="$BACKGROUND_ANALYSIS"
        [[ -z "$force_run" ]] && force_run="false"

        if [[ -z "$output_dir" ]]; then
            if [[ -n "$GLOBAL_OUTPUT_DIR" ]]; then
                output_dir="${GLOBAL_OUTPUT_DIR}/${task_id}"
            else
                output_dir="${input_folder}/output/${task_id}"
            fi
        fi

        if [[ -z "$input_folder" || -z "$chromosome" || -z "$center" ]]; then
            echo "${task_id},${input_folder},${chromosome},${center},${window},${step},${background},${output_dir},${force_run},FAILED,missing required fields" >> "$batch_summary"
            fail_count=$((fail_count + 1))
            continue
        fi

        echo "[BATCH] Running ${task_id}: ${input_folder} ${chromosome}:${center}"

        if run_single_task "$input_folder" "$chromosome" "$center" "$window" "$step" "$background" "$output_dir" "$force_run"; then
            echo "${task_id},${input_folder},${chromosome},${center},${window},${step},${background},${output_dir},${force_run},OK,completed" >> "$batch_summary"
            success_count=$((success_count + 1))
        else
            echo "${task_id},${input_folder},${chromosome},${center},${window},${step},${background},${output_dir},${force_run},FAILED,see pipeline.log" >> "$batch_summary"
            fail_count=$((fail_count + 1))
        fi
    done < "$batch_file"

    echo "[BATCH] Finished. Success=${success_count}, Failed=${fail_count}"
    echo "[BATCH] Summary: $batch_summary"
}

if [[ -n "$BATCH_CSV" ]]; then
    run_batch "$BATCH_CSV"
else
    if [[ -z "$INPUT_FOLDER" || -z "$CHROMOSOME" || -z "$CENTER_POSITION" ]]; then
        echo "[ERROR] Missing required arguments for single task mode."
        usage
    fi

    if [[ -z "$OUTPUT_DIR" ]]; then OUTPUT_DIR="${INPUT_FOLDER}/output"; fi
    run_single_task "$INPUT_FOLDER" "$CHROMOSOME" "$CENTER_POSITION" "$WINDOW_SIZE" "$STEP_SIZE" "$BACKGROUND_ANALYSIS" "$OUTPUT_DIR" "$FORCE_RUN"
fi