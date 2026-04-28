import requests

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def resolve_voc(name: str) -> dict:
    # Pas 1: obtenir el CID
    url_cid = f"{PUBCHEM_URL}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    r = requests.get(url_cid, timeout=10)
    r.raise_for_status()
    cid = r.json()["IdentifierList"]["CID"][0]
    print(f"  CID trobat: {cid}")

    # Pas 2: obtenir propietats
    url_props = f"{PUBCHEM_URL}/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES,InChIKey/JSON"
    r2 = requests.get(url_props, timeout=10)
    r2.raise_for_status()
    props = r2.json()["PropertyTable"]["Properties"][0]

    # PubChem pot retornar el camp amb noms diferents segons el compost
    smiles = (props.get("IsomericSMILES")
              or props.get("CanonicalSMILES")
              or props.get("SMILES")
              or props.get("ConnectivitySMILES"))

    inchikey = props.get("InChIKey")

    print(f"  InChIKey: {inchikey}")
    print(f"  SMILES: {smiles}")

    return {
        "name": name,
        "smiles": smiles,
        "inchikey": inchikey
    }

if __name__ == "__main__":
    result = resolve_voc("1-octen-3-ol")
    print(result)