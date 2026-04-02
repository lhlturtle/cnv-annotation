# CNV注释系统

根据染色体位置输出cytoband注释和遗传咨询术语。

**参考基因组**: hg19

## 数据来源

- **CytoBand**: UCSC hg19 cytoBand
- **基因注释**: NCBI Gene Info (193,000+ 基因)
- **临床区域**: OMIM, ClinGen, DECIPHER

## 功能

- 输入：染色体号、起始位置、终止位置、Gain/Loss、ratio值
- 输出1：染色体cytoband结果（ISCN表示）
- 输出2：遗传咨询术语（包含基因注释和临床意义）

## 示例

输入：`chr11, 81500000, 135000000, Loss, 0.521`

输出1：`del(11)(q14.1-q25)(53.5Mb)`

输出2：该样本检测到11号染色体q14.1-q25区域存在约53.5Mb的缺失...

## 快速启动

```bash
cd cnv-annotation
docker-compose up -d
```

访问：http://localhost:8502

## 技术栈

- Python 3.11
- Streamlit
- Docker
