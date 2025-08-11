FROM python:3.11-slim
WORKDIR /app
RUN apt-get update && apt-get install -y \
    build-essential \
    libpq-dev \
    curl \
    libxrender-dev \
    libxext-dev \
    && rm -rf /var/lib/apt/lists/*
RUN pip install mcpo uv 
RUN pip install mcp 

COPY . /app
RUN pip install .

ENTRYPOINT ["/bin/bash", "-c"," mcpo --host 0.0.0.0 --port 8000 --  mcp-server-template"]