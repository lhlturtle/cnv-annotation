# CNV 注释项目

## 项目概述
- **类型**：遗传分析工具
- **功能**：输入染色体位置，输出cytoband和遗传咨询术语
- **输入**：染色体号、起始位置、终止位置、Gain/Loss、ratio值
- **输出**：
  1. 染色体cytoband结果（如 del(11)(q14.1-q25)(53.5Mb)）
  2. 遗传咨询术语

## 技术方案
- Python + Streamlit (Web界面)
- 本地cytoband数据库
- 基因注释数据库（OMIM, ClinGen, DECIPHER）
- Docker部署

## 数据来源
- CytoBand: UCSC Genome Browser
- 基因注释: 本地数据库

## 进度
- [x] 项目创建
- [ ] cytoband数据
- [ ] 基因注释
- [ ] 界面开发
- [ ] 测试
