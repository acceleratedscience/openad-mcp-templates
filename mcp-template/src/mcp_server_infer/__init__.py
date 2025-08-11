from . import server
import asyncio


def main():
    asyncio.run(server.server.run_stdio_async())
