#!/usr/bin/env python
# coding: utf-8

"""
Author: Khoi Van 2020

This script takes a CAS number and look up its pKa (dissociation constant) 
from Pubchem and return it if found; return None, otherwise

Change notes:
- 2020-02-14:
- Instead of just returning the result, return the result combo as json dict

"""

import re
import sys
import traceback
import xml.etree.ElementTree as ET
from typing import Optional

import pandas as pd
import pubchempy as pcp  # https://pubchempy.readthedocs.io/en/latest/guide/gettingstarted.html
import requests

from .classify import classify


debug = False


def pka_lookup_pubchem(identifier, namespace=None, domain='compound') -> Optional[str]:
    """
    Lookup the pKa value of a compound from PubChem.

    This function uses the PubChem API to find the pKa value of a compound based on a given identifier.
    The identifier can be a name, CAS number, SMILES, InChI, or InChIKey. The function also supports debugging mode.

    Args:
        identifier (str): The identifier for the compound. This can be a chemical name, CAS number, SMILES, InChI, or InChIKey.
        namespace (str, optional): The namespace of the identifier. This can be 'name', 'cas', 'smiles', 'inchi', or 'inchikey'.
                                   If not provided, the function will attempt to classify the identifier.
        domain (str, optional): The domain to search within. Default is 'compound'.

    Returns:
        Optional[str]: A dictionary with the pKa value and other relevant information if found, otherwise None.

    Raises:
        ValueError: If an exact match for the identifier is not found on PubChem.
        RuntimeError: If the compound or pKa value is not found in PubChem.
    """
    global debug

    if len(sys.argv) == 2 and sys.argv[1] in ['--debug=True', '--debug=true', '--debug', '-d']:
        debug = True

    # Identify lookup source (Pubchem in this case)
    lookup_source = 'Pubchem'

    try:
        headers = {
            'user-agent': 'Mozilla/5.0 (X11; CentOS; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/73.0.3683.75 Safari/537.36'}

        cids = []
        identifier_type = ''

        if not namespace:
            identifier_type = classify(identifier)

            if identifier_type in ['smiles', 'inchi', 'inchikey']:
                lookup = pcp.get_cids(identifier, namespace=identifier_type)
                if lookup:
                    cids.append(lookup[0])
            else:
                lookup = pcp.get_cids(identifier, namespace='name')
                if lookup:
                    cids.append(lookup[0])
        elif namespace == 'cas':
            cids = pcp.get_cids(identifier, namespace='name')
        else:
            cids = pcp.get_cids(identifier, namespace=namespace)

        if not cids:
            lookup = pcp.get_cids(identifier, namespace='name')
            if lookup:
                cids.append(lookup[0])
            identifier_type = namespace

        if len(cids) > 0:
            cid = cids[0]
            exact_match = True
            synonyms = pcp.get_synonyms(cid)[0]['Synonym'] or []

            returned_cas = ''
            for synonym in synonyms:
                cas_nr = re.search(r'^\d{2,7}-\d{2}-\d$', synonym)
                if cas_nr:
                    cas_nr = cas_nr.group()
                    returned_cas = cas_nr
                    break

            lookup_result = pcp.get_properties(['inchi', 'inchikey',
                                                'canonical_smiles', 'isomeric_smiles',
                                                'iupac_name'],
                                               cid)

            if identifier_type == 'cas':
                exact_match = identifier in synonyms
            elif identifier_type in ['inchi', 'inchikey']:
                if identifier_type == 'inchi':
                    exact_match = (identifier == lookup_result[0].get('InChI', False))
                elif identifier_type == 'inchikey':
                    exact_match = (identifier == lookup_result[0].get('InChIKey', False))

            if not exact_match:
                if debug:
                    print(f'Exact match between input and Pubchem return value? {identifier in synonyms}')
                raise ValueError('This is not an exact match on Pubchem!')

            pka_lookup_result_xml = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/XML?heading=Dissociation+Constants'.format(cid)

            r = requests.get(pka_lookup_result_xml, headers=headers, timeout=15)
            if r.status_code == 200 and len(r.history) == 0:
                tree = ET.fromstring(r.text)
                info_node = tree.find('.//*{http://pubchem.ncbi.nlm.nih.gov/pug_view}Information')
                original_source = info_node.find('{http://pubchem.ncbi.nlm.nih.gov/pug_view}Reference').text
                pka_result = info_node.find('.//*{http://pubchem.ncbi.nlm.nih.gov/pug_view}String').text
                pka_result = re.sub(r'^pKa = ', '', pka_result)

                core_result = {
                    'source': lookup_source,
                    'Pubchem_CID': str(cid),
                    'pKa': pka_result,
                    'reference': original_source,
                    'Substance_CASRN': returned_cas,
                }
                extra_info = lookup_result[0]
                extra_info.pop('CID', None)

                result = {**core_result, **extra_info}
                s = pd.Series(result)
                s = s.rename({
                    'CanonicalSMILES': 'Canonical_SMILES',
                    'IsomericSMILES': 'Isomeric_SMILES',
                    'IUPACName': 'IUPAC_Name'
                })
                result = s.to_dict()
                return result

            else:
                raise RuntimeError('pKa not found in Pubchem.')

        else:
            raise RuntimeError('Compound not found in Pubchem.')

    except Exception as error:
        if debug:
            traceback_str = ''.join(traceback.format_exception(etype=type(error), value=error, tb=error.__traceback__))
            print(traceback_str)

        return None


if __name__ == "__main__":
    import pprint as pp
    # from pprint import pprint as print
    print = pp.PrettyPrinter(indent=1, width=80).pprint

    cas_nr = '64-19-7'    # acetic acid   >>> pKa = 4.76 at 25 °C
    # cas_nr = '75-75-2'    # methanesulfonic acid   >>> pKa = -1.86
    # cas_nr = '2950-43-8'    # Hydroxylamine-O-sulfonic acid, no result
    # cas_nr = '2687-12-9'
    # print(pka_lookup_pubchem(cas_nr))
    print(pka_lookup_pubchem(cas_nr, 'cas'))



    # smiles_string = 'C1=CC(=CC=C1F)S'
    # smiles_string = 'OC=1(N(N=C(C=1)C2(=CC=CC=C2))C)'
    smiles_string = 'OC1=CC=CC=C1'

    # # Look up pKa using pka_lookup_pubchem():
    # print(f'pKa from Pubchem using smiles: {pka_lookup_pubchem(smiles_string)}')
    # print('pKa from Pubchem using smiles:')
    # print(pka_lookup_pubchem(smiles_string, "smiles"))

    inchi_string = 'InChI=1S/C6H6S/c7-6-4-2-1-3-5-6/h1-5,7H'
    # inchi_string = 'InChI=1S/C10H10N2O/c1-12-10(13)7-9(11-12)8-5-3-2-4-6-8/h2-7,13H,1H3'    # this is NOT exact match from Pubchem return search
    # inchi_string = 'InChI=1S/C10H10N2O/c1-12-10(13)7-9(11-12)8-5-3-2-4-6-8/h2-7,11H,1H3'    # this is an exact match from Pubchem return search

    # print(f'pKa from Pubchem using smiles: {pka_lookup_pubchem(inchi_string)}')
    # print(f'pKa from Pubchem using smiles: {pka_lookup_pubchem(inchi_string, "smiles")}')    # this function call has wrong 'namespace' (should be 'inchi', not 'smiles'). Therefore, return Pubchem CID even though it is not an exact match.
    # print(f'pKa from Pubchem using smiles:')
    # print(pka_lookup_pubchem(inchi_string, "inchi"))

    inchikey_string = 'OKKJLVBELUTLKV-UHFFFAOYSA-N'    # methanol
    # print(f'pKa from Pubchem using InChIKey:')
    # print(pka_lookup_pubchem(inchikey_string, "inchikey"))