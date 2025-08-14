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
import threading

REQUESTOR = OpenadAPI()
REQUESTOR.request("?")  # initialising handle
REQUESTOR_2 = OpenadAPI()
REQUESTOR_2.request("?")
REQUESTOR_3 = OpenadAPI()
REQUESTOR_3.request("?")

# from .process_creds import place_models
import json
import requests

server = FastMCP("SMI-TED-PFAS")

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



async def get_alternative_names_a_molecule(molecule: str) -> dict:
    """
    get alternative names or identifiers for a molecules.
        accepts a molecule inchi string, smiles string or name as input
        keys in the dictionary:
               - 'identifiers' Dictionary of the known identifiers smiles, inchi and molecular formula
               - 'synonynms' list of other names it is known as e.g. drug names or chemical names


        for example:
                     {
    "identifiers": {
      "name": "ethanol",
      "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
      "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
      "canonical_smiles": "CCO",
      "isomeric_smiles": "CCO",
      "molecular_formula": "C2H6O",
      "cid": 702
    },
    "synonyms":["alcohol","hand sanitiser"] }
    """
    global REQUESTOR
    mol = REQUESTOR.request(f"export mol {molecule}")
    row = {}
    row["identifiers"] = {}
    for key in mol["identifiers"]:
        if key == "record":
            continue
        if mol["identifiers"][key]:
            row["identifiers"][key] = mol["identifiers"][key]

    row["synonyms"] = mol["synonyms"]
    dump = json.dumps(row, indent=2)
    return row


def possible_smiles(smiles: str) -> bool:
    """
    Verify is a string *could* be a SMILES definition.

    Parameters
    ----------
    smiles: str
        The SMILES string to validate
    """
    return bool(re.search(r"[BCNOFPSI](?:[a-df-z0-9#=@+%$:\[\]\(\)\\\/\.\-])*", smiles, flags=re.I))


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
async def get_smi_properties(smiles: list, properties: list):
    """get the one or more of the following properties:
    [ qm8-e1-cam, qm8-e1-cc2, qm8-e1-pbe0, qm8-e2-cam, qm8-e2-cc2, qm8-e2-pbe0, qm8-f1-cam, qm8-f1-cc2, qm8-f1-pbe0, qm8-f2-cam, qm8-f2-cc2, qm8-f2-pbe0, qm9-alpha, qm9-cv, qm9-g298, qm9-gap, qm9-h298, qm9-lumo, qm9-homo, qm9-mu, qm9-r2, qm9-u0, qm9-u298, qm9-zpve,bace, bbbp, biodegradability, clintox, esol, freesolv, hiv, iupac, ld50, lipo, logkow, sider, tox21 ]

    These results are based on the SMILES Ted Foundation model

    Note: The word score or result in a request may refer to the result of the analysis.
    convert all identifed properties to lowecase

    Args:
        smiles (list): list of smiles strings e.g. ["CCO", "CC(O)OC"]
        properties (list) : list of properties from the following list in uppercase from the following [ qm8-e1-cam, qm8-e1-cc2, qm8-e1-pbe0, qm8-e2-cam, qm8-e2-cc2, qm8-e2-pbe0, qm8-f1-cam, qm8-f1-cc2, qm8-f1-pbe0, qm8-f2-cam, qm8-f2-cc2, qm8-f2-pbe0, qm9-alpha, qm9-cv, qm9-g298, qm9-gap, qm9-h298, qm9-lumo, qm9-homo, qm9-mu, qm9-r2, qm9-u0, qm9-u298, qm9-zpve,bace, bbbp, biodegradability, clintox, esol, freesolv, hiv, iupac, ld50, lipo, logkow, sider, tox21 ]
    "

        Returns:
        dictionary DataFrame containing target smiles, properties and property results

     Notes:
    If calling this from a UI like Open UI Ensure when calling a tool that any list parameter is passed as a list e.g ["Apple","Orange"] Not a string containing a list  " [\"Apple\",\'Orange'\" "

    Here are descriptions of the supported properties
    QM9 Properties Summary:
| Property | Description |
|----------|-------------|
| **qm9-alpha** | Isotropic polarizability (unit: Bohr^3) |
| **qm9-cv** | Heat capacity at constant volume at 298.15K (unit: cal/(mol*K)) |
| **qm9-g298** | Gibbs free energy at 298.15K (unit: Hartree) |
| **qm9-gap** | Gap between HOMO and LUMO (unit: Hartree) |
| **qm9-h298** | Enthalpy at 298.15K (unit: Hartree) |
| **qm9-lumo** | Lowest unoccupied molecular orbital energy (unit: Hartree) |
| **qm9-homo** | Highest occupied molecular orbital energy (unit: Hartree) |
| **qm9-mu** | Dipole moment (unit: Debye) |
| **qm9-r2** | Electronic spatial extent (unit: Bohr^2) |
| **qm9-u0** | Internal energy at 0K (unit: Hartree) |
| **qm9-u298** | Internal energy at 298.15K (unit: Hartree) |
| **qm9-zpve** | Zero point vibrational energy (unit: Hartree) |

### MoleculeNet Properties Summary:
| Property | Description |
|----------|-------------|
| **bace** | Inhibition of human beta secretase 1 |
| **bbbp** | Blood brain barrier penetration |
| **biodegradability** | Predicting the biodegradability of compounds |
| **clintox** | Toxicity data of FDA-approved drugs and those that fail clinical trials |
| **esol** | Delaney water solubility data for organics |
| **freesolv** | Hydration free energy |
| **hiv** | Inhibition of HIV viral replication |
| **iupac** | IUPAC (description incomplete) |
| **ld50** | Toxicity modeling (description incomplete) |
| **lipo** | Octonol/water distribution coefficient |
| **logkow** | (description incomplete) |
| **sider** | Side Effect Resource. Market drugs and their adverse drug reactions/side effects |
| **tox21** | Toxicity measurements on 12 different targets |

### QM8 summaries Summary:
Transition Energies (E)
These parameters represent the energy required for an electron to transition from the ground state (S0) to an excited singlet state (S1 or S2).

qm8-e1-cam: S0 -> S1 transition energy calculated with LR-TDCAM-B3LYP/def2TZVP.

qm8-e1-cc2: S0 -> S1 transition energy calculated with RI-CC2/def2TZVP.

qm8-e1-pbe0: S0 -> S1 transition energy calculated with LR-TDPBE0/def2??VP. 

qm8-e2-cam: S0 -> S2 transition energy calculated with LR-TDCAM-B3LYP/def2TZV. (This might be a slight variation from def2TZVP.)

qm8-e2-cc2: S0 -> S2 transition energy calculated with RI-CC2/def2TZVP.

qm8-e2-pbe0: S0 -> S2 transition energy calculated with LR-TDPBE0/def2??VP. 

Oscillator Strengths (f)
These parameters indicate the probability of an electron transitioning from the ground state (S0) to an excited singlet state (S1 or S2). Higher values suggest a stronger transition.

qm8-f1-cam: S0 -> S1 oscillator strength calculated with LR-TDCAM-B3LYP/def2TZV. (Note: The original description had a minor discrepancy, but this refers to the S1 state.)

qm8-f1-cc2: S0 -> S1 oscillator strength calculated with RI-CC2/def2TZVP.

qm8-f1-pbe0: S0 -> S1 oscillator strength calculated with LR-TDPBE0/def2??VP. 

qm8-f2-cam: S0 -> S2 oscillator strength calculated with LR-TDCAM-B3LYP/def2TZV.

qm8-f2-cc2: S0 -> S2 oscillator strength calculated with RI-CC2/def2TZVP.

qm8-f2-pbe0: S0 -> S2 oscillator strength calculated with LR-TDPBE0/def2??VP. 
    """


    global REQUESTOR
    # requestor = OpenadAPI()
    qm9 = [
        "qm9-alpha",
        "qm9-cv",
        "qm9-g298",
        "qm9-gap",
        "qm9-h298",
        "qm9-lumo",
        "qm9-homo",
        "qm9-mu",
        "qm9-r2",
        "qm9-u0",
        "qm9-u298",
        "qm9-zpve",
    ]
    qm8 = [
        "qm8-e1-cam",
        "qm8-e1-cc2",
        "qm8-e1-pbe0",
        "qm8-e2-cam",
        "qm8-e2-cc2",
        "qm8-e2-pbe0",
        "qm8-f1-cam",
        "qm8-f1-cc2",
        "qm8-f1-pbe0",
        "qm8-f2-cam",
        "qm8-f2-cc2",
        "qm8-f2-pbe0",
    ]
    molnet = [
        "bace",
        "bbbp",
        "biodegradability",
        "clintox",
        "esol",
        "freesolv",
        "hiv",
        "iupac",
        "ld50",
        "lipo",
        "logkow",
        "sider",
        "tox21",
    ]
    qm9_target = []
    qm8_target = []
    molnet_target = []
    for aprop in properties:
        prop = aprop.lower()
        if prop in qm9:
            qm9_target.append(prop)
        elif prop in qm8:
            qm8_target.append(prop)
        elif prop in molnet:
            molnet_target.append(prop)

    logger.info("for Smiles :" + str(smiles) + "getting the properties :" + str(properties))
    try:
        valid_smiles_list = [valid_smiles(smile) for smile in smiles]
        smiles_list = []
        for mol in valid_smiles_list:
            if mol is not None:
                smiles_list.append(mol)
        logger.info(f"bmfm_sm get molecule property {properties} for {smiles_list} ")
        if len(smiles_list) == 0:
            return "unable to verify list of smiles"

        # result = REQUESTOR.request(f"bmfm_sm get molecule property {properties} for {smiles_list} ")
        threads = []
        results = []
        if len(qm8_target) > 0:
            t1 = threading.Thread(
                target=lambda: results.append(
                    REQUESTOR.request(
                        f"smi get molecule property {qm8_target} for {smiles_list} ",
                    )
                )
            )
            threads.append(t1)
        if len(qm9_target) > 0:
            t2 = threading.Thread(
                target=lambda: results.append(
                    REQUESTOR_3.request(
                        f"smi get molecule property {qm9_target} for {smiles_list} ",
                    )
                )
            )
            threads.append(t2)
        if len(molnet_target) > 0:
            t3 = threading.Thread(
                target=lambda: results.append(
                    REQUESTOR_2.request(f"smi get molecule property {molnet_target} for {smiles_list} ")
                )
            )
            threads.append(t3)

        for x in threads:
            x.start()

        for x in threads:
            x.join()
        logger.info("results\n " + str(results))
        if len(results) == 0:
            return "no positive result"
        else:
            result = pd.concat(results, ignore_index=True)
            return result.to_dict()
    except Exception as e:
        logger.info("unable to return a result see below for details")
        logger.error(str(e))
        logger.error(str(traceback.format_exc().splitlines()))

        return "no positive result"


@server.tool()
async def get_pfas_classification(smiles: list):
    """
    Function: determine is a molecule is classified as a PFAS or PFOS forever chemical under given statuatory regulations
    classifications  OECD 2021 , ECHA 2021 struct and EPA 2023

    Args:
    smiles (list): of valid smiles strings e.g. ["C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O", "CCO"]

    Returns:
    dictionary with classifications


    Notes:
    If calling this from a UI like Open UI Ensure when calling a tool that any list parameter is passed as a list e.g ["Apple","Orange"] Not a string containing a list  " [\"Apple\",\'Orange'\" "

    """
    global REQUESTOR
    # requestor = OpenadAPI()

    try:
        valid_smiles_list = [valid_smiles(smile) for smile in smiles]
        smiles_list = []
        for mol in valid_smiles_list:
            if mol is not None:
                smiles_list.append(mol)
        if len(smiles_list) == 0:
            return "unable to verify list of smiles"
        logger.info(f"pfas get molecule property [ OECD_2021 , ECHA_2021_struct, EPA_2023 ]  for {smiles_list} ")

        result = REQUESTOR.request(
            f"pfas get molecule property [ OECD_2021 , ECHA_2021_struct, EPA_2023 ]  for {smiles_list} "
        )
        logger.info(result)
        return result.to_dict()
    except Exception as e:
        logger.info("unable to return a result see below for details")
        logger.error(str(e))
        logger.error(str(traceback.format_exc().splitlines()))

        return "no positive result"


async def main():
    asyncio.run(server.run_stdio_async())


if __name__ == "__main__":

    main()
