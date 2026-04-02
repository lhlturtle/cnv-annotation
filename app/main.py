"""
CNV注释系统 - Web界面
"""
import streamlit as st
import sys
sys.path.insert(0, '/app')

from annotation import annotate_cnv

# 页面配置
st.set_page_config(
    page_title="CNV注释系统",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 CNV注释系统")
st.markdown("输入染色体位置，获取cytoband注释和遗传咨询术语")

# 输入表单
with st.expander("📥 CNV输入", expanded=True):
    col1, col2 = st.columns([1, 1])
    
    with col1:
        chr_num = st.text_input("染色体号", value="11", help="如: 11 或 chr11")
        start_pos = st.number_input("起始位置(bp)", value=81500000, min_value=0, step=1000000)
        end_pos = st.number_input("终止位置(bp)", value=135000000, min_value=0, step=1000000)
    
    with col2:
        cn_type = st.selectbox("变异类型", ["Loss", "Gain"], index=0)
        ratio = st.number_input("Ratio值", value=0.521, min_value=0.0, max_value=2.0, format="%.6f")
    
    if st.button("🔍 注释", type="primary"):
        if chr_num and start_pos < end_pos:
            try:
                iscn_result, counseling_result = annotate_cnv(
                    chr_num, start_pos, end_pos, cn_type, ratio
                )
                
                # 显示结果
                st.success("注释完成！")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("### 📍 ISCN表示")
                    st.code(iscn_result, language="text")
                    st.caption("染色体核型表示法")
                
                with col2:
                    st.markdown("### 📋 遗传咨询术语")
                    st.text_area("遗传咨询报告", value=counseling_result, height=300)
                
                # 复制按钮
                st.button("📋 复制结果")
                
            except Exception as e:
                st.error(f"注释失败: {e}")
        else:
            st.warning("请检查输入参数")

# 示例
st.markdown("---")
st.subheader("📝 使用示例")

example_data = [
    {"chr": "11", "start": 81500000, "end": 135000000, "type": "Loss", "ratio": 0.521},
]

for ex in example_data:
    if st.button(f"示例: chr{ex['chr']} {ex['start']}-{ex['end']} {ex['type']}"):
        try:
            iscn_result, counseling_result = annotate_cnv(
                ex['chr'], ex['start'], ex['end'], ex['type'], ex['ratio']
            )
            st.session_state.example_iscn = iscn_result
            st.session_state.example_counseling = counseling_result
        except Exception as e:
            st.error(f"处理失败: {e}")

if 'example_iscn' in st.session_state:
    st.markdown("### 📍 示例结果")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**ISCN表示**")
        st.code(st.session_state.example_iscn)
    with col2:
        st.markdown("**遗传咨询术语**")
        st.text(st.session_state.example_counseling)

# 说明
st.markdown("---")
with st.expander("ℹ️ 使用说明"):
    st.markdown("""
    ### 输入参数
    - **染色体号**: 1-22, X, Y
    - **起始/终止位置**: 基因组坐标(bp)
    - **变异类型**: Loss(缺失) 或 Gain(重复)
    - **Ratio值**: 拷贝数分析的比值
    
    ### 输出说明
    1. **ISCN表示**: 国际标准染色体命名法
    2. **遗传咨询术语**: 包含基因注释和临床意义解读
    
    ### 数据来源
    - CytoBand: UCSC Genome Browser
    - 基因注释: OMIM, ClinGen, DECIPHER
    """)
