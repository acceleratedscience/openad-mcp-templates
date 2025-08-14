from mcp.server.fastmcp import FastMCP
import logging
import asyncio
from pydantic import BaseModel
import sys
import pandas as pd
from openad import OpenadAPI
from rdkit import Chem
import re
import traceback

REQUESTOR = OpenadAPI()
# from .process_creds import place_models
import json
import requests

server = FastMCP("BMFM-PFAS")

from dotenv import load_dotenv, dotenv_values

# loading variables from .env file
load_dotenv()
# place_models()


handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s")
handler.setFormatter(formatter)
# Create a logger
logger = logging.getLogger("Demo-model-server")
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


def valid_smiles(smiles: str):
    """
    Verify if a string is valid SMILES definition.

    Parameters
    ----------
    smiles: str
        The SMILES string to validate
    """
    logger.info("checking validity of smiles string: " + str(smiles))
    if not smiles or not possible_smiles(smiles):
        logger.error("Smiles not possible:" + str(smiles))
        return None

    try:

        m = Chem.MolFromSmiles(smiles, sanitize=False)  # pylint: disable=no-member
    except Exception:  # pylint: disable=broad-exception-caught
        logger.error("Error converting smiles to mol:" + str(smiles))
        return None

    if m is None:
        return None
    else:
        try:
            Chem.SanitizeMol(m)  # pylint: disable=no-member
        except Exception:  # pylint: disable=broad-exception-caught
            logger.error("Unable to sanitize mol:" + str(smiles))
            return None
    try:
        result = Chem.MolToSmiles(m, canonical=True)
    except:
        logger.error("Unable to convert mol to smiles :" + str(smiles))
        return None
    return result


@server.tool()
async def get_bmfm_properties(smiles: list, properties: list):
    """get the one or more of the following properties [BACE, BBBP, CLINTOX, ESOL, FREESOLV, HIV, LIPOPHILICITY, MUV, QM7, SIDER, TOX21, TOXCAST] for a list of Smiles strings

    the word score or result in a request may refer to the result of the analysis.

    Args:
        smiles (list): list of smiles strings e.g. ["CCO", "CC(O)OC"]
        properties (list) : list of properties from the following list in uppercase from the following ["BACE", "BBBP", "CLINTOX", "ESOL", "FREESOLV", "HIV", "LIPOPHILICITY", "MUV", "QM7", "SIDER", "TOX21", "TOXCAST"]
    "

        Returns:
        dictionary DataFrame containing target smiles, properties and property results

     Notes:
    If calling this from a UI like Open UI Ensure when calling a tool that any list parameter is passed as a list e.g ["Apple","Orange"] Not a string containing a list  " [\"Apple\",\'Orange'\" "
    """

    global REQUESTOR
    # requestor = OpenadAPI()
    valid_properties = [
        "BACE",
        "BBBP",
        "CLINTOX",
        "ESOL",
        "FREESOLV",
        "HIV",
        "LIPOPHILICITY",
        "MUV",
        "QM7",
        "SIDER",
        "TOX21",
        "TOXCAST",
    ]
    target_properties = []
    for aprop in properties:
        prop = aprop.upper()
        if prop in valid_properties:
            target_properties.append(prop)

    logger.info("for Smiles :" + str(smiles) + "getting the properties :" + str(target_properties))
    try:
        valid_smiles_list = [valid_smiles(smile) for smile in smiles]
        smiles_list = []
        for mol in valid_smiles_list:
            if mol is not None:
                smiles_list.append(mol)
        if len(target_properties) == 0 or len(smiles_list) == 0:
            return "no positive result, incorrect smiles or property lists"

        logger.info(f"bmfm_sm get molecule property {target_properties} for {smiles_list} ")
        if len(smiles_list) == 0:
            return "unable to verify list of smiles"
        result = REQUESTOR.request(f"bmfm_sm get molecule property {target_properties} for {smiles_list} ")
        logger.info(str(result))
        return result.to_dict()
    except Exception as e:
        logger.info("unable to return a result see below for details")
        logger.error(str(e))
        logger.error(str(traceback.format_exc().splitlines()))

        return "no positive result"


@server.tool()
async def get_protein_solubility(protein_strings: list):
    """
    Function: calculate the solubility (sol) for a list of protein strings

    the protein string is provided in FASTA format and the drug molecule ins smiles string format.

    Args:
    protein_strings (list): Protein FASTA Strings  e.g. ["NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC","NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC","NLMKCTESTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"]

    Returns:
    dictionary DataFrame containing target proteins and solubility score

    Notes:
    If calling this from a UI like Open UI Ensure when calling a tool that any list parameter is passed as a list e.g ["Apple","Orange"] Not a string containing a list  " [\"Apple\",\'Orange'\" "

    """
    global REQUESTOR
    # requestor = OpenadAPI()
    try:
        result = REQUESTOR.request(f"bmfm_pm get protein property sol for {protein_strings} ")
        logger.info((f"bmfm_pm get protein property sol for {protein_strings} "))
        logger.info(str(result))
        return result.to_dict()
    except:

        return "no positive result"


@server.tool()
async def get_dti_interaction(protein_strings: list, drug_smiles: str):
    """
    Function: calculate the Drug to target (dti)  binding score for a list of proteins and drug molecules

    calculate the drug target interaction (dti) score between a list of protein strings and a drug smiles string

    the protein string is provided in FASTA format and the drug molecule ins smiles string format.

    Args:
    protein_strings (list): Protein FASTA Strings  e.g. ["NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC","NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC","NLMKCTESTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"]

    drug_smiles (str): Drug smiles string e.g. 'CCO'


    Returns:
    dictionary DataFrame containing target protein and binding affinity score

    Notes:
    If calling this from a UI like Open UI Ensure when calling a tool that any list parameter is passed as a list e.g ["Apple","Orange"] Not a string containing a list  " [\"Apple\",\'Orange'\" "
    """
    global REQUESTOR

    # requestor = OpenadAPI()
    try:
        drug_smile = valid_smiles(drug_smiles)
        logger.info(f"bmfm_pm get protein property dti for {protein_strings} using ( drug_smiles='{drug_smile}' )")
        result = REQUESTOR.request(
            f"bmfm_pm get protein property dti for {protein_strings} using ( drug_smiles='{drug_smile}' )"
        )
        logger.info(str(result))
        return result.to_dict()
    except:

        return "no positive result"


@server.tool()
async def search_fasta_sequence(fasta_string: str):
    """
    Search the RCSB PDB for a given Protein FASTA string and returns its details

    Parameters:
        fasta_string: str


    Returns:

    """
    global REQUESTOR
    sequence_type = "protein"
    return_first = True

    # https://search.rcsb.org/#search-example-3
    search_dict = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.9,
                "sequence_type": sequence_type,
                "value": fasta_string,
            },
        },
        "request_options": {"scoring_strategy": "sequence"},
        "return_type": "entry",
    }

    # Make search_dict URL safe
    search_str = json.dumps(search_dict)
    search_str = search_str.replace(" ", "")
    search_str = encode_uri_component(search_str)

    # Search for the sequence in the PDB
    search_url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={search_str}"
    search_response = requests.get(search_url)
    search_results = search_response.json() if search_response.status_code == 200 else {}

    # Error handling
    if not search_results.get("result_set"):
        return False, "No matching PDB entries found."

    # Return the first result if it's a 100% match
    if len(search_results["result_set"]) > 0 and search_results["result_set"][0].get("score") == 1:
        pdb_id = search_results["result_set"][0]["identifier"]
        return fetch_pdb_file(pdb_id)
    else:
        return False, search_results["result_set"]


def encode_uri_component(string):
    from urllib.parse import quote

    return quote(string.encode("utf-8"), safe="~()*!.'")


def fetch_pdb_file(pdb_id, file_format="cif"):
    """
    Fetch a PDB file by its ID from rscb.org.

    Parameters:
        pdb_id: str
            The PDB ID to fetch.
        file_format: 'cif' | 'pdb'
            The format of the PDB file to fetch.

    Returns:
        success: bool
            Whether the request was successful.
        file_data: str
            The content of the (cif) file.
    """

    # Fetch the file data.
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
    pdb_response = requests.get(pdb_url)

    if pdb_response.status_code != 200:
        return False, "Failed to retrieve PDB data."

    file_data = pdb_response.text
    return True, file_data


async def main():
    asyncio.run(server.run_stdio_async())


if __name__ == "__main__":

    main()
