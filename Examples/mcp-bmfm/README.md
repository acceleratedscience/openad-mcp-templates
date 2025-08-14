# This is a simple template

### Note
 for this to run you will need to add an Key to access the OpenAD Gateway. If you do not have a key you can start the respective bmfm services and catalog them in the dockerfile
<br>
to run the template using the Echo function:

### 1: run `podman build -t mcp-server-template:latest .`

### 2: run `podman compose up`

### Note:
To apply your own code to the template simply add your required packages to the `pyproject.toml` 
<br>
Then edit the `src/mcp_server_infer/server.py` and add your tool functions with the decorator `@server.tool`<br>

This service will run on port 8000 the OpenAI standard MCP protocol supported by Webui etc, it will also support on 8001 streamable_http.

recreate the container image and r-run compose... your services will be available on `127.0.0.1:8000/openapi.json`

