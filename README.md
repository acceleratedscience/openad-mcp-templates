
# Simple Template MCP
=======================

## How It Works
---------------

The stack uses the MCP defined from FASTAPI and MCPO Openapi.

### Creating Your MCP
---------------------

1. Under the template, you can find the `src/mcp_server_infer/server.py` library.
2. This is where you place your Library code and define your MCPs. The standard MCP definition is used as per [https://github.com/modelcontextprotocol/python-sdk](https://github.com/modelcontextprotocol/python-sdk).
3. We have chosen to use the main MCP libraries as they were found to be the most reliable.

### Running Your MCP Server
---------------------------

1. Simply install your project and run your MCP server by running:
    ```
    python -m mcp_server_infer
    ```
2. By default, the MCP service will run on port 8000 using stdio.
3. If you wish to use another protocol, you can use the following command options:
    * For SSE protocol: `python -m mcp_server_infer sse`
    * For Streamable HTTP protocol: `python -m mcp_server_infer shttp`
    * For stdio protocol: `python -m mcp_server_infer stdio`

### Serving as OpenAPI
------------------------
We are using MCPO from https://github.com/open-webui/mcpo


1. To serve this in addition to the MCP server, invoke the MCPO proxy service:
    ```
    mcpo --port 8000 --host 0.0.0.0 --server-type "streamable_http" -- http://127.0.0.1:8001/mcp
    ```
2. Alternatively, you can run stdio and port 8000 with OpenAPI format by running:
    ```
    mcpo mcp_server_infer
    ```



## Directories

Both the Examples and template provide compose files and Docker builds. the docker file uses the given entrypoint.sh whch starts first the shttp protocol server then the mcpo proxy.<br>

### mcp-template
This Directory has a simple Echo examples of an MCP qeith ***Openapi*** interface exposed on port `8000` and ***streamble http*** on `8001`


### Examples

This directory shows how to create inferences running against the OpenAD Biomedical Model Protein and Small Molecule Inferences.

Current setup requires access to the OpenAD , gateway. Subsequent iteractions will provide single Instance desktop deployments using podman