#!/bin/sh

set -e

HOST=0.0.0.0
PORT=8001
TIMEOUT=1200
START=$(date +%s)

echo "Starting mcp-server-infer service on port $PORT..."

# add the openad credentials to the mcp-server-infer service
python /app/src/mcp_server_infer/process_creds.py
# Start the mcp-server-infer service in the background
python -m mcp_server_infer shttp &

MCP_PID=$!

echo "Waiting for mcp-server-infer service at $HOST:$PORT to become available..."

# Wait for the service to be ready
while ! nc -z $HOST $PORT; do
  sleep 1
  NOW=$(date +%s)
  if [ $((NOW - START)) -ge $TIMEOUT ]; then
    echo "Timeout after $TIMEOUT seconds waiting for $HOST:$PORT"
    kill $MCP_PID 2>/dev/null
    exit 1
  fi
done

echo "mcp-server-infer service is ready at $HOST:$PORT. Starting mcpo proxy..."

# Start the mcpo proxy service the blow starts it a "streamable_http" this can be changed to sse to support sse
mcpo --port 8000 --host 0.0.0.0 --server-type "streamable_http" -- http://127.0.0.1:8001/mcp
