import server
import asyncio


if __name__ == "__main__":
    asyncio.run(server.get_bmfm_properties(["CCO", "CC(O)OC"], ["BACE", "ESOL"]))

    asyncio.run(
        server.get_dti_interaction(
            ["NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"],
            "C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O",
        )
    )

    asyncio.run(server.get_protein_solubility(["NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"]))

    asyncio.run(server.get_dti_interaction(["NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"], "CC(O)OC"))

    asyncio.run(server.get_pfas_classification(["C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O"]))
    asyncio.run(server.get_pfas_classification(["C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O"]))

    asyncio.run(server.get_pfas_classification(["C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O", "CCO"]))
