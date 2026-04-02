FROM python:3.11-slim

WORKDIR /app

# 安装依赖
COPY app/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 复制代码
COPY app/ /app/

# 创建数据目录
RUN mkdir -p /app/data

# 暴露端口
EXPOSE 8502

# 启动命令
CMD ["streamlit", "run", "main.py", "--server.address=0.0.0.0", "--server.port=8502"]
