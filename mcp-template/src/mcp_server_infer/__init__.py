import os
import sys
from . import server
import asyncio


def main():
    if len(sys.argv) > 1:
        server_type = sys.argv[
            1
        ].lower()  # Convert to lowercase for consistent checking
    else:
        # Default to 'stdio' or 'default' if no argument is provided
        server_type = "stdio"

    if server_type.lower() == "sse":
        server.server.settings.port = 8001
        server.server.settings.host = "0.0.0.0"
        asyncio.run(server.server.run_sse_async())
    elif server_type.lower() == "shttp":
        server.server.settings.port = 8001
        server.server.settings.host = "0.0.0.0"
        asyncio.run(server.server.run_streamable_http_async())
    else:
        # server.server.settings.port = 8002
        asyncio.run(server.server.run_stdio_async())
