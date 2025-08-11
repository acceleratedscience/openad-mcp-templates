from mcp.server.fastmcp import FastMCP
import logging
import asyncio
from pydantic import BaseModel

import re
import traceback


server = FastMCP("model-server-template")


handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s")
handler.setFormatter(formatter)
# Create a logger
logger = logging.getLogger("Demo-model-server")
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


@server.tool()
async def mcp_echo(echo_message: str):
    """
    Echos the result

    Args:
    echo_message (str): The message to be echoed back

    Returns:
    String
    """

    return "the ECHO ...." + echo_message


async def main():
    asyncio.run(server.run_stdio_async())


if __name__ == "__main__":

    main()
