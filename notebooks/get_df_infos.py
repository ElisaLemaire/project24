import pandas as pd
from pubchemprops.pubchemprops import get_cid_by_name, get_first_layer_props, get_second_layer_props
import urllib.error
import urllib.parse
from pka_lookup import pka_lookup_pubchem
import pubchempy as pcp
import re
import json

#Finds the pKa using the code of Khoi Van.
def find_pka(inchikey_string):
    text_pka = pka_lookup_pubchem(inchikey_string, "inchikey")
    if text_pka is not None and 'pKa' in text_pka:
        pKa_value = text_pka['pKa']
        return pKa_value
    else:
        return None

def find_boiling_point(name):
    text_dict = get_second_layer_props(str(name), ['Boiling Point', 'Vapor Pressure'])
    
    Boiling_point_values = []
    #finds all celsius
    pattern_celsius = r'([-+]?\d*\.\d+|\d+) °C'
    pattern_F = r'([-+]?\d*\.\d+|\d+) °F'
    
    for item in text_dict['Boiling Point']:
        # Check if the item has a key 'Value' and 'StringWithMarkup'
        if 'Value' in item and 'StringWithMarkup' in item['Value']:
            # Access the 'String' key inside the nested dictionary
            string_value = item['Value']['StringWithMarkup'][0]['String']
            match_celsius = re.search(pattern_celsius, string_value)
            if match_celsius:
                celsius = float(match_celsius.group(1))
                Boiling_point_values.append(celsius)

            #Search for Farenheit values, if found: converts farenheit to celsius before adding to the list
            match_F = re.search(pattern_F, string_value)
            if match_F:
                fahrenheit_temp = float(match_F.group(1))
                celsius_from_F = round(((fahrenheit_temp - 32) * (5/9)), 2)
                Boiling_point_values.append(celsius_from_F)

    #get the mean value
    Boiling_temp = round((sum(Boiling_point_values) / len(Boiling_point_values)), 2)
    return Boiling_temp

"""
This code takes as input a list of compound written like: acetone, water. The code allows spaces, wrong names and unknown pubchem names.
Then it iterates through each of them to find if they exist on pubchem, and if they do,
then 'CID', 'MolecularFormula', 'MolecularWeight', 'InChI', 'InChIKey', 'IUPACName', 'XLogP', 'pKa',  and 'BoilingPoint' is added into a list and then a data frame.
The code takes time as find_pka(inchikey_string) and find_boiling_point(name) request URL to find the string on the Pubchem page, then extract it using regex. 
The boiling point is a mean of all the values (references) found.
"""

def get_df_properties():
    compound_name = input("Enter the name of the compound like this: water, ethanol, methanol =")
    compound_list = [compound.strip() for compound in compound_name.split(',')]
    
    compound_properties = []
    valid_properties = []
    for compound_name in compound_list:
        compound_name_encoded = urllib.parse.quote(compound_name.strip())
        try: 
            first_data = get_first_layer_props(compound_name_encoded, ['MolecularFormula', 'MolecularWeight', 'InChI', 'InChIKey', 'IUPACName', 'XLogP'])
            compound_info = {}
            for prop in ['CID', 'MolecularFormula', 'MolecularWeight', 'InChI', 'InChIKey', 'IUPACName', 'XLogP']:
                compound_info[prop] = first_data.get(prop)
            #print(first_data)
            
            #adds pKa value
            pka_value = find_pka(first_data['InChIKey'])
            if pka_value is not None:
                compound_info['pKa'] = pka_value
            else:
                pass
            
            #adds boiling point, solubility
            compound_info['BoilingPoint'] = find_boiling_point(compound_name_encoded)
    
            # When every property has been added to compound_info, add to the properties. This makes sure all properties have the right keys
            compound_properties.append(compound_info)
        
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f'{compound_name} not found on PubChem')
            else:
                print(f'An error occurred: {e}')
    
    for prop in compound_properties:
        if isinstance(prop, dict):
            valid_properties.append(prop)
    df = pd.DataFrame(valid_properties)
    # Set the property names from the first dictionary as column headers
    if len(valid_properties) > 0:
        df = df.reindex(columns=valid_properties[0].keys())

    #print(df)
    return(df)