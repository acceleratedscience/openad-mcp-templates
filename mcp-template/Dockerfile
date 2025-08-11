FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    libpq-dev \
    curl \
    libxrender-dev \
    libxext-dev \
    busybox \
    netcat-openbsd \
    && rm -rf /var/lib/apt/lists/*

RUN pip install mcpo uv 
#RUN pip install mcp openad

COPY . .
RUN pip install .

RUN chmod +x /app/entrypoint.sh

EXPOSE 8000
EXPOSE 8001

ENTRYPOINT [ "/app/entrypoint.sh" ]
