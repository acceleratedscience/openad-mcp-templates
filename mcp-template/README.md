# This is a simple template
<br>
to run the template using the Echo function:

### 1: run `podman build -t mcp-server-template:latest .`

### 2: run `podman compose up`

### Note:
To apply your own code to the template simply add your required packages to the `pyproject.toml` 
<br>
Then edit the `src/mcp_server_infer/server.py` and add your tool functions with the decorator `@server.tool`

recreate the container image and r-run compose... your services will be available on `127.0.0.1:8000/openapi.json`

To access them from Open Webui simply add the above url to your Seettings

![Screenshot 2025-05-16 at 4 05 20 PM](https://github.ibm.com/Accelerated-Discovery/ad-mcp-templates/assets/225313/4099a0f5-d1b0-458f-b7b1-2618d60d0928)
