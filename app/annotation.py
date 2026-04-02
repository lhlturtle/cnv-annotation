"""
CNV注释模块 - 染色体位置到cytoband和遗传咨询术语转换
支持 hg19 基因组版本
"""
import json
import re
from pathlib import Path

# 数据路径
DATA_DIR = Path(__file__).parent.parent / "data"

def load_cytoband_data():
    """加载cytoband数据 - hg19"""
    # 优先使用hg19
    cytoband_file = DATA_DIR / "cytoband_hg19.txt"
    data = {}
    
    if not cytoband_file.exists():
        cytoband_file = DATA_DIR / "cytoband_hg38.txt"
    
    if not cytoband_file.exists():
        cytoband_file = DATA_DIR / "cytoband.txt"
    
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
                stain = parts[4] if len(parts) > 4 else ''
                
                if chr_name not in data:
                    data[chr_name] = []
                data[chr_name].append({
                    'start': start,
                    'end': end,
                    'band': band,
                    'stain': stain
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
    """简化cytoband名称，如 p15.1 -> (p, 15.1)"""
    match = re.match(r'([pq])(\d+\.?\d*)', band)
    if match:
        return match.group(1), match.group(2)
    return band, ''

def get_arm(chr_num, position):
    """根据位置判断染色体臂"""
    # 对于大多数染色体，p臂在60Mb之前
    # 但这对某些染色体不准确，这里用简化判断
    p_arm_end = {
        '1': 125000000, '2': 93300000, '3': 91000000, '4': 50600000,
        '5': 48400000, '6': 60800000, '7': 59900000, '8': 45600000,
        '9': 43000000, '10': 40200000, '11': 53700000, '12': 35800000,
        '13': 17600000, '14': 17600000, '15': 19000000, '16': 36600000,
        '17': 25100000, '18': 18700000, '19': 26500000, '20': 27500000,
        '21': 13200000, '22': 14700000, 'X': 58600000, 'Y': 10400000
    }
    threshold = p_arm_end.get(str(chr_num), 60000000)
    return 'p' if position < threshold else 'q'

def format_iscn_notation(chr_num, start_pos, end_pos, cn_type):
    """生成ISCN格式的CNV表示"""
    start_band, end_band = get_cytoband_range(chr_num, start_pos, end_pos)
    size_mb = (end_pos - start_pos) / 1_000_000
    
    # 确定CNV类型
    if cn_type == 'Loss':
        prefix = 'del'
    elif cn_type == 'Gain':
        prefix = 'dup'
    else:
        prefix = 'cnv'
    
    if not start_band or not end_band:
        # 如果找不到精确的cytoband，使用坐标
        start_arm = get_arm(chr_num, start_pos)
        end_arm = get_arm(chr_num, end_pos)
        if start_arm == end_arm:
            return f"{prefix}({chr_num})({start_arm}{int(start_pos/1_000_000)}-{end_arm}{int(end_pos/1_000_000)})({size_mb:.1f}Mb)"
        else:
            return f"{prefix}({chr_num})(p{int(start_pos/1_000_000)}-q{int(end_pos/1_000_000)})({size_mb:.1f}Mb)"
    
    # 简化band名称
    start_arm, start_sub = simplify_cytoband(start_band)
    end_arm, end_sub = simplify_cytoband(end_band)
    
    # 生成ISCN表示
    if start_arm == end_arm:
        if start_sub == end_sub:
            return f"{prefix}({chr_num})({start_band})({size_mb:.1f}Mb)"
        else:
            return f"{prefix}({chr_num})({start_band}-{end_band})({size_mb:.1f}Mb)"
    else:
        return f"{prefix}({chr_num})({start_band}-{end_band})({size_mb:.1f}Mb)"

def find_matching_region(chr_num, start_pos, end_pos, annotation_data):
    """查找匹配的基因注释区域"""
    chr_num = str(chr_num).replace('chr', '')
    
    for region_key, region_data in annotation_data.items():
        if region_key.startswith('_'):
            continue
        if region_data.get('chromosome', '').replace('chr', '') == chr_num:
            region_start = region_data.get('start', 0)
            region_end = region_data.get('end', 0)
            # 检查是否有重叠
            if start_pos <= region_end and end_pos >= region_start:
                return region_key, region_data
    return None, None

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
    
    # 查找匹配的注释区域
    region_key, region_data = find_matching_region(chr_num, start_pos, end_pos, annotation_data)
    
    # 构建遗传咨询文本
    report = []
    
    # 第一段：基础检测结果
    report.append(f"该样本检测到{chr_num}号染色体")
    
    if start_band and end_band:
        report.append(f"{start_band}-{end_band}")
    elif start_band:
        report.append(f"{start_band}")
    
    report.append(f"（chr{chr_num}:{start_pos}-{end_pos}）区域存在约{size_mb:.1f}Mb的{cnv_type_cn}，")
    report.append(f"拷贝数为{copy_num}。")
    
    # 第二段：基因概述
    if region_key and region_data:
        omim_genes = region_data.get('omim_genes', {})
        clingen_genes = region_data.get('clingen_genes', {})
        
        if omim_genes:
            report.append(f"\n受检样本该{cnv_type_cn}区域包含大量重要的可编码蛋白质基因。")
            report.append(f"\n其中OMIM疾病相关基因有：")
            
            gene_list = []
            for gene, info in list(omim_genes.items())[:8]:
                gene_str = f"{gene}"
                if 'omim_id' in info:
                    gene_str += f" ({info['omim_id']})"
                if 'description' in info and info['description']:
                    gene_str += f"（{info['description']}）"
                gene_list.append(gene_str)
            
            report.append("，".join(gene_list))
            report.append("等。")
        
        # 第三段：ClinGen证据
        if clingen_genes:
            report.append(f"\n\n该{cnv_type_cn}区域在ClinGen数据库中包含了已知的单倍体剂量敏感性基因/区域，")
            report.append("其临床意义明确：")
            
            for gene, info in clingen_genes.items():
                if 'omim_id' in info:
                    report.append(f"\n{gene}（{info['omim_id']}）基因：")
                else:
                    report.append(f"\n{gene}基因：")
                
                if 'disease' in info:
                    report.append(f"与{info['disease']}相关。")
                
                if 'description' in info:
                    report.append(f"{info['description']}。")
                
                if 'phenotype' in info:
                    report.append(f"其主要临床表型包括：{info['phenotype']}。")
                
                if 'haploinsufficiency_score' in info:
                    report.append(f"\n  单倍体剂量评分(HI-Score): {info['haploinsufficiency_score']}")
        
        # 第四段：DECIPHER案例
        decipher_cases = region_data.get('decipher_cases', [])
        if decipher_cases:
            report.append("\n\n疾病数据库DECIPHER中有与受检样本该缺失区域存在重叠的记录，表型具有明确关联。")
            report.append("\n其中不乏明确致病或可能致病案例：")
            
            for case in decipher_cases[:3]:
                report.append(f"\nPatient: {case.get('patient_id', 'Unknown')}（{case.get('inheritance', '未知')}）：")
                report.append(f"表现为{case.get('phenotype', '未知表型')}。")
        
        # 第五段：综合评估
        path_terms = {
            'Pathogenic': '致病性变异（Pathogenic, P）',
            'Likely Pathogenic': '可能致病性变异（Likely Pathogenic, LP）',
            'VUS': '临床意义不明变异（Variant of Uncertain Significance, VUS）',
            'Benign': '可能良性变异（Likely Benign, LB）',
            'Likely Benign': '可能良性变异（Likely Benign, LB）'
        }
        
        path_score = region_data.get('pathogenicity', 'VUS')
        acmg_class = region_data.get('acmg_classification', 'VUS')
        
        report.append(f"\n\n综合该CNV的变异类型、位置、包含基因功能以及数据库情况，")
        
        if start_band and end_band:
            report.append(f"受检样本该{start_band}-{end_band}缺失为")
        else:
            report.append(f"受检样本该区域缺失为")
        
        report.append(f"{path_terms.get(path_score, path_score)}。")
        
        if acmg_class and acmg_class != 'VUS':
            report.append(f"\n参考ACMG分类：{acmg_class}")
    
    else:
        # 没有匹配的注释区域
        report.append(f"\n\n该区域涉及多个基因，功能分析需要进一步研究。")
        report.append(f"\n建议：结合家族史、临床表现及其他检测结果进行综合分析。")
        report.append(f"\n可参考DECIPHER、ClinGen等数据库进一步查询该区域的临床意义。")
    
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
    print("=" * 60)
    print("测试用例：chr11, 81500000, 135000000, Loss, 0.521")
    print("=" * 60)
    result = annotate_cnv("chr11", 81500000, 135000000, "Loss", 0.521)
    print("\nISCN表示:", result[0])
    print("\n遗传咨询术语:")
    print(result[1])
