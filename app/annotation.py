"""
CNV注释模块 - 染色体位置到cytoband和遗传咨询术语转换
"""
import json
import re
from pathlib import Path

# 数据路径
DATA_DIR = Path(__file__).parent.parent / "data"

def load_cytoband_data():
    """加载cytoband数据"""
    cytoband_file = DATA_DIR / "cytoband.txt"
    data = {}
    
    with open(cytoband_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chr_name = parts[0].replace('chr', '')
                start = int(parts[1])
                end = int(parts[2])
                band = parts[3]
                
                if chr_name not in data:
                    data[chr_name] = []
                data[chr_name].append({
                    'start': start,
                    'end': end,
                    'band': band,
                    'stain': parts[4] if len(parts) > 4 else ''
                })
    
    return data

def load_gene_annotation():
    """加载基因注释数据"""
    annotation_file = DATA_DIR / "gene_annotation.json"
    with open(annotation_file, 'r', encoding='utf-8') as f:
        return json.load(f)

def get_cytoband_for_position(chr_name, position):
    """获取指定位置对应的cytoband"""
    cytoband_data = load_cytoband_data()
    chr_name = str(chr_name).replace('chr', '')
    
    if chr_name not in cytoband_data:
        return None
    
    for band_info in cytoband_data[chr_name]:
        if band_info['start'] <= position <= band_info['end']:
            return band_info['band']
    
    return None

def get_cytoband_range(chr_num, start_pos, end_pos):
    """获取一段区域的cytoband范围"""
    cytoband_data = load_cytoband_data()
    chr_num = str(chr_num).replace('chr', '')
    
    if chr_num not in cytoband_data:
        return None, None
    
    start_band = None
    end_band = None
    
    for band_info in cytoband_data[chr_num]:
        if band_info['start'] <= start_pos <= band_info['end']:
            start_band = band_info['band']
        if band_info['start'] <= end_pos <= band_info['end']:
            end_band = band_info['band']
    
    return start_band, end_band

def simplify_cytoband(band):
    """简化cytoband名称，如 p15.1 -> p15.1"""
    # 移除多余的点数
    match = re.match(r'([pq])(\d+\.?\d*)', band)
    if match:
        return match.group(1), match.group(2)
    return band, ''

def format_iscn_notation(chr_num, start_pos, end_pos, cn_type):
    """生成ISCN格式的CNV表示"""
    start_band, end_band = get_cytoband_range(chr_num, start_pos, end_pos)
    
    if not start_band or not end_band:
        # 如果找不到精确的cytoband，使用坐标
        size_mb = (end_pos - start_pos) / 1_000_000
        arm = 'q' if start_pos > 60_000_000 else 'p'  # 简化的臂判断
        return f"{'del' if cn_type == 'Loss' else 'dup'}({chr_num})({arm}{int(start_pos/1_000_000)}-{arm}{int(end_pos/1_000_000)})({size_mb:.1f}Mb)"
    
    # 简化band名称
    _, start_sub = simplify_cytoband(start_band)
    _, end_sub = simplify_cytoband(end_band)
    
    start_arm = 'q' if 'q' in start_band else 'p'
    end_arm = 'q' if 'q' in end_band else 'p'
    
    # 计算大小
    size_mb = (end_pos - start_pos) / 1_000_000
    
    # 确定CNV类型
    if cn_type == 'Loss':
        prefix = 'del'
    elif cn_type == 'Gain':
        prefix = 'dup'
    else:
        prefix = 'cnv'
    
    # 生成ISCN表示
    if start_arm == end_arm:
        if start_sub == end_sub:
            return f"{prefix}({chr_num})({start_band})({size_mb:.1f}Mb)"
        else:
            return f"{prefix}({chr_num})({start_arm}{start_sub}-{end_sub})({size_mb:.1f}Mb)"
    else:
        return f"{prefix}({chr_num})({start_band}-{end_band})({size_mb:.1f}Mb)"

def generate_genetic_counseling(chr_num, start_pos, end_pos, cn_type, ratio):
    """生成遗传咨询术语"""
    
    # 获取cytoband信息
    start_band, end_band = get_cytoband_range(chr_num, start_pos, end_pos)
    size_mb = (end_pos - start_pos) / 1_000_000
    
    # 确定CNV类型中文
    if cn_type == 'Loss':
        cnv_type_cn = "缺失"
        copy_num = 1
    else:
        cnv_type_cn = "重复"
        copy_num = 3
    
    # 加载基因注释
    annotation_data = load_gene_annotation()
    
    # 查找匹配的注释
    region_key = None
    if start_band and end_band:
        # 构造key，如 "11q14.1-q25"
        arm = 'q' if 'q' in start_band else 'p'
        region_key = f"{chr_num}{arm}{start_band.split(arm)[1] if arm in start_band else start_band}"
        # 简化搜索
        for key in annotation_data:
            if arm in key and start_band.replace(arm, '') in key:
                region_key = key
                break
    
    # 构建遗传咨询文本
    report = []
    report.append(f"该样本检测到{chr_num}号染色体")
    
    if start_band and end_band:
        report.append(f"{start_band}-{end_band}")
    
    report.append(f"（chr{chr_num}:{start_pos}-{end_pos}）区域")
    report.append(f"存在约{size_mb:.1f}Mb的{cnv_type_cn}，")
    report.append(f"拷贝数为{copy_num}。")
    report.append(f"\n受检样本该{cnv_type_cn}区域包含大量重要的可编码蛋白质基因。")
    
    # 添加基因信息
    if region_key and region_key in annotation_data:
        anno = annotation_data[region_key]
        
        if 'omim_genes' in anno and anno['omim_genes']:
            omim_genes = anno['omim_genes']
            gene_list = [f"{gene} ({info['omim_id']})" for gene, info in omim_genes.items()]
            report.append(f"\n其中OMIM疾病相关基因有：")
            report.append(", ".join(gene_list[:5]))
            report.append("等。")
        
        if 'clingen_genes' in anno and anno['clingen_genes']:
            report.append(f"\n该{cnv_type_cn}区域在ClinGen数据库中包含了已知的单倍体剂量敏感性基因/区域，")
            report.append("其临床意义明确：")
            
            for gene, info in anno['clingen_genes'].items():
                report.append(f"\n{gene}（{gene}）区域：包含{gene}（{info['omim_id']}）基因。")
                if 'note' in info:
                    report.append(f"该基因的单倍体剂量不足与{info['description']}相关。")
    
    # 添加DECIPHER案例
    if region_key and region_key in annotation_data:
        anno = annotation_data[region_key]
        if 'decipher_cases' in anno:
            report.append("\n\n疾病数据库DECIPHER中有与受检样本该缺失区域存在重叠的记录，表型具有明确关联。")
            report.append("\n其中不乏明确致病或可能致病案例：")
            
            for case in anno['decipher_cases']:
                report.append(f"\nPatient: {case['patient_id']}（{case['inheritance']}）：表现为{case['phenotype']}。")
    
    # 综合评估
    if region_key and region_key in annotation_data:
        anno = annotation_data[region_key]
        if 'pathogenicity' in anno:
            path_terms = {
                'Pathogenic': '致病性变异（Pathogenic, P）',
                'Likely Pathogenic': '可能致病性变异（Likely Pathogenic, LP）',
                'VUS': '临床意义不明变异（VUS）',
                'Benign': '可能良性变异（Likely Benign, LB）'
            }
            report.append(f"\n\n综合该CNV的变异类型、位置、包含基因功能以及数据库情况，")
            report.append(f"受检样本该{start_band if start_band else chr_num}区域{cnv_type_cn}为{path_terms.get(anno['pathogenicity'], anno['pathogenicity'])}。")
    
    return "".join(report)

def annotate_cnv(chr_num, start_pos, end_pos, cn_type, ratio):
    """
    CNV注释主函数
    
    参数:
        chr_num: 染色体号 (如 "11" 或 "chr11")
        start_pos: 起始位置 (bp)
        end_pos: 终止位置 (bp)
        cn_type: CNV类型 ("Loss" 或 "Gain")
        ratio: 拷贝数比值
    
    返回:
        tuple: (ISCN表示, 遗传咨询术语)
    """
    chr_num = str(chr_num).replace('chr', '')
    
    # 输出1: ISCN表示
    iscn = format_iscn_notation(chr_num, int(start_pos), int(end_pos), cn_type)
    
    # 输出2: 遗传咨询术语
    counseling = generate_genetic_counseling(chr_num, int(start_pos), int(end_pos), cn_type, ratio)
    
    return iscn, counseling

if __name__ == "__main__":
    # 测试
    result = annotate_cnv("chr11", 81500000, 135000000, "Loss", 0.5212835047122986)
    print("ISCN表示:", result[0])
    print("\n遗传咨询术语:")
    print(result[1])
