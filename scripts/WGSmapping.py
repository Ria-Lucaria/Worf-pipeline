import subprocess
import csv
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
import os
import sys
import argparse
from datetime import datetime

# ==========================================
# 工具函数
# ==========================================

def check_and_fix_chrom_name(bam_file, target_chrom):
    """检查 BAM 文件头，自动修正染色体名称"""
    print(f"[DEBUG] 正在检查染色体名称匹配: 输入 '{target_chrom}' vs BAM文件...")
    try:
        cmd = ['samtools', 'view', '-H', bam_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        bam_chroms = set()
        for line in result.stdout.splitlines():
            if line.startswith('@SQ'):
                parts = line.split('\t')
                for p in parts:
                    if p.startswith('SN:'):
                        bam_chroms.add(p[3:])
        
        if target_chrom in bam_chroms:
            return target_chrom
        if f"chr{target_chrom}" in bam_chroms:
            print(f"[WARN] 自动修正染色体名称: {target_chrom} -> chr{target_chrom}")
            return f"chr{target_chrom}"
        if target_chrom.startswith('chr'):
            no_chr = target_chrom.replace('chr', '')
            if no_chr in bam_chroms:
                print(f"[WARN] 自动修正染色体名称: {target_chrom} -> {no_chr}")
                return no_chr
                
        print(f"[ERROR] BAM文件中未找到染色体: {target_chrom} (且无法自动修正)")
        return target_chrom 
    except Exception as e:
        print(f"[WARN] 无法读取 BAM Header: {e}")
        return target_chrom

def get_bam_chrom_length(bam_file, chrom):
    """获取指定染色体的长度"""
    try:
        cmd = ['samtools', 'view', '-H', bam_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.split('\n'):
            if line.startswith('@SQ') and f'SN:{chrom}' in line:
                for part in line.split('\t'):
                    if part.startswith('LN:'):
                        return int(part[3:])
    except Exception:
        pass
    return None

def get_depth_histogram(bam_file, chrom, start, end):
    """使用 samtools depth 快速获取逐碱基覆盖度"""
    positions = []
    depths = []
    region = f"{chrom}:{start}-{end}"
    cmd = ['samtools', 'depth', '-r', region, bam_file]
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            cols = line.strip().split('\t')
            if len(cols) < 3: continue
            positions.append(int(cols[1]))
            depths.append(int(cols[2]))
        process.wait()
    except Exception as e:
        print(f"[ERROR] samtools depth 失败: {e}")
    return positions, depths

def get_read_intervals(bam_file, chrom, start, end):
    """获取单条Reads的区间坐标"""
    intervals = []
    region = f"{chrom}:{start}-{end}"
    cmd = ['samtools', 'view', '-F', '4', bam_file, region]
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            cols = line.split('\t')
            if len(cols) < 10: continue
            
            read_start = int(cols[3])
            read_len = len(cols[9]) 
            read_end = read_start + read_len
            
            if read_end < start or read_start > end:
                continue
            intervals.append((read_start, read_end))
        process.wait()
    except Exception as e:
        print(f"[ERROR] 读取 Reads 失败: {e}")
    return intervals

def get_hit_reads(bam_file, chrom, start, end):
    """提取窗口内命中 reads 的详细信息"""
    hit_reads = []
    region = f"{chrom}:{start}-{end}"
    cmd = ['samtools', 'view', '-F', '4', bam_file, region]

    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in process.stdout:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 11:
                continue

            read_name = cols[0]
            flag = int(cols[1])
            read_chrom = cols[2]
            read_start = int(cols[3])
            mapq = int(cols[4])
            cigar = cols[5]
            seq = cols[9]

            read_len = len(seq)
            read_end = read_start + read_len

            if read_end < start or read_start > end:
                continue

            overlap_start = max(start, read_start)
            overlap_end = min(end, read_end)
            overlap_bp = max(0, overlap_end - overlap_start + 1)

            strand = '-' if (flag & 16) else '+'

            hit_reads.append({
                'read_name': read_name,
                'chromosome': read_chrom,
                'read_start': read_start,
                'read_end': read_end,
                'read_length': read_len,
                'mapq': mapq,
                'flag': flag,
                'strand': strand,
                'cigar': cigar,
                'overlap_start': overlap_start,
                'overlap_end': overlap_end,
                'overlap_bp': overlap_bp,
            })
        process.wait()
    except Exception as e:
        print(f"[ERROR] 读取命中 Reads 详情失败: {e}")

    return hit_reads

def write_hit_reads_csv(hit_reads, out_csv):
    fieldnames = [
        'read_name', 'chromosome', 'read_start', 'read_end', 'read_length',
        'mapq', 'flag', 'strand', 'cigar', 'overlap_start', 'overlap_end', 'overlap_bp'
    ]
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in hit_reads:
            writer.writerow(row)
    print(f"✅ 生成命中 Reads 明细: {out_csv}")

def write_summary(summary_path, sample_name, chrom, center, start_pos, end_pos, depth_vals, hit_reads):
    total_hits = len(hit_reads)
    unique_reads = len({x['read_name'] for x in hit_reads})
    mean_depth = float(np.mean(depth_vals)) if len(depth_vals) > 0 else 0.0
    max_depth = int(np.max(depth_vals)) if len(depth_vals) > 0 else 0
    nonzero_depth_bp = int(np.sum(np.array(depth_vals) > 0)) if len(depth_vals) > 0 else 0
    window_size_bp = end_pos - start_pos + 1
    coverage_ratio = (nonzero_depth_bp / window_size_bp) if window_size_bp > 0 else 0.0

    with open(summary_path, 'w') as f:
        f.write(f"sample_name: {sample_name}\n")
        f.write(f"region: {chrom}:{start_pos}-{end_pos}\n")
        f.write(f"center: {center}\n")
        f.write(f"window_size_bp: {window_size_bp}\n")
        f.write(f"total_hit_reads: {total_hits}\n")
        f.write(f"unique_hit_reads: {unique_reads}\n")
        f.write(f"mean_depth: {mean_depth:.4f}\n")
        f.write(f"max_depth: {max_depth}\n")
        f.write(f"covered_bp: {nonzero_depth_bp}\n")
        f.write(f"coverage_ratio: {coverage_ratio:.6f}\n")

    print(f"✅ 生成任务总结: {summary_path}")

def greedy_stacking(intervals, gap=1):
    intervals.sort(key=lambda x: x[0])
    rows = [] 
    packed_data = [] 
    
    for start, end in intervals:
        placed = False
        length = end - start
        
        for i in range(len(rows)):
            if rows[i] + gap < start:
                rows[i] = end 
                packed_data.append((start, length, i))
                placed = True
                break
        
        if not placed:
            rows.append(end)
            packed_data.append((start, length, len(rows) - 1))
            
    return packed_data, len(rows)

# ==========================================
# 绘图函数 (核心修改部分)
# ==========================================

def plot_pileup(packed_data, total_rows, chrom, start, end, title, filename, center_pos=None):
    """
    绘制线性堆叠图
    Update 1: 使用 tab20 色板进行分层着色
    Update 2: 强制 Y 轴最小高度，使 Reads 看起来更扁平
    """
    # 保持画布大小不变
    fig, ax = plt.subplots(figsize=(14, 6)) 
    
    # 1. 设置颜色
    # 找到 plot_pileup 函数中的这一行
    # cmap = plt.get_cmap('tab20')  <-- 删除这行
    
    # 换成这行：
    cmap = plt.get_cmap('Set1')

    patches = []
    face_colors = []
    
    for (s, l, r) in packed_data:
        # Rectangle((x, y), width, height)
        rect = Rectangle((s, r), l, 0.8) # 高度保持0.8，留0.2空隙
        patches.append(rect)
        # 根据行号 r 循环取色
        # face_colors.append(cmap(r % 20))
        face_colors.append(cmap(r % 9))  # Set1 色板有9种颜色
    
    # 2. 创建集合并上色
    coll = PatchCollection(patches, edgecolor='none', alpha=1.0)
    coll.set_facecolor(face_colors) # 应用颜色列表
    ax.add_collection(coll)
    
    # 3. 设置坐标轴 (实现"压扁"效果的关键)
    ax.set_xlim(start, end)
    
    # 逻辑：如果层数很少(比如5层)，强制把Y轴拉到50，这样5层占的空间就很小，柱子就扁了
    # 如果层数很多(比如100层)，就用实际层数+缓冲
    Y_MIN_VISUAL_LIMIT = 50 
    y_limit = max(total_rows + 2, Y_MIN_VISUAL_LIMIT)
    ax.set_ylim(0, y_limit)
    
    # 装饰
    ax.set_xlabel(f"Position on {chrom} (bp)", fontsize=12)
    ax.set_ylabel("Stacked Depth (Layer)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # 标记中心点
    if center_pos:
        ax.axvline(x=center_pos, color='#333333', linestyle='--', linewidth=1, alpha=0.8)
        # 标签也相应往上移
        ax.text(center_pos, y_limit * 1.01, f" {center_pos:,}", color='#333333', 
                transform=ax.transData, fontweight='bold', fontsize=10)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"✅ 生成堆叠图: {filename}")

def plot_local_histogram(pos_list, depth_list, chrom, start, end, title, filename):
    """绘制局部覆盖度柱状图"""
    if not pos_list: return
    plt.figure(figsize=(14, 4))
    plt.fill_between(pos_list, depth_list, color='gray', alpha=0.6, step="mid")
    plt.xlim(start, end)
    plt.xlabel(f"Position on {chrom} (bp)")
    plt.ylabel("Depth")
    plt.title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"✅ 生成对照组柱状图: {filename}")

# ==========================================
# 主程序
# ==========================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True)
    parser.add_argument('--chromosome', required=True)
    parser.add_argument('--center', type=int, required=True)
    parser.add_argument('--step', type=int, default=100000)
    parser.add_argument('--window', type=int, default=10000)
    parser.add_argument('--background', default='false')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    bam_file = args.bam
    raw_chrom = args.chromosome
    center = args.center
    half_window = args.window
    out_dir = args.output
    
    os.makedirs(out_dir, exist_ok=True)
    basename = os.path.splitext(os.path.basename(bam_file))[0]
    sample_name = basename.replace('_aligned.sorted', '').replace('_aligned', '')

    chrom = check_and_fix_chrom_name(bam_file, raw_chrom)
    
    start_pos = max(1, center - half_window)
    end_pos = center + half_window
    print(f"[INFO] 分析窗口: {chrom}:{start_pos:,}-{end_pos:,}")

    # 统一时间戳与命名前缀：包含 -c, -p, -w 与时间戳
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # 清理可能的路径不安全字符（理论上 chrom 是安全的，这里保持一致性）
    chrom_safe = str(chrom).replace('/', '_')
    prefix = f"{sample_name}_c{chrom_safe}_p{center}_w{half_window}_{timestamp}"

    # Step 1: Histogram Check
    depth_pos, depth_vals = get_depth_histogram(bam_file, chrom, start_pos, end_pos)
    if len(depth_vals) > 0:
        hist_file = os.path.join(out_dir, f"{prefix}_hist.png")
        plot_local_histogram(depth_pos, depth_vals, chrom, start_pos, end_pos, 
                            f"Coverage Check: {chrom}:{center}", hist_file)

    # Step 2: 命中 reads 统计输出
    hit_reads = get_hit_reads(bam_file, chrom, start_pos, end_pos)
    hit_csv = os.path.join(out_dir, f"{prefix}_hit_reads.csv")
    summary_file = os.path.join(out_dir, f"{prefix}_summary.txt")
    write_hit_reads_csv(hit_reads, hit_csv)
    write_summary(summary_file, sample_name, chrom, center, start_pos, end_pos, depth_vals, hit_reads)

    # Step 3: Pile-up Plot
    intervals = [(x['read_start'], x['read_end']) for x in hit_reads]
    
    if len(intervals) > 0:
        if len(intervals) > 30000:
            print("[WARN] Reads > 30000, subsampling...")
            intervals = intervals[::2]

        packed, max_rows = greedy_stacking(intervals, gap=1)
        print(f"    堆叠完成，最大层数: {max_rows}")
        
        out_png = os.path.join(out_dir, f"{prefix}_pileup.png")
        plot_pileup(packed, max_rows, chrom, start_pos, end_pos,
                   f"Reads Pile-up: {chrom}:{center}",
                   out_png, center_pos=center)
    else:
        print("[WARN] 没有 Reads，跳过绘图。")

if __name__ == "__main__":
    main()