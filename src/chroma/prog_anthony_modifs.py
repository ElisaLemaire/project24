import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re
import pandas as pd
from pubchemprops.pubchemprops import get_cid_by_name, get_first_layer_props, get_second_layer_props
import urllib.error
import urllib.parse
from pka_lookup import pka_lookup_pubchem
import re
import json

# Function to add a molecule to the mixture
def add_molecule():
    element = mixture_entry.get()
    mixture.append(element)
    mixture_listbox.insert(tk.END, element.strip())
    return mixture, mixture_listbox


"""
This code takes as input a list of compound written like: acetone, water. The code allows spaces, wrong names and unknown pubchem names.
Then it iterates through each of them to find if they exist on pubchem, and if they do,
then 'CID', 'MolecularFormula', 'MolecularWeight', 'InChI', 'InChIKey', 'IUPACName', 'XLogP', 'pKa',  and 'BoilingPoint' is added into a list and then a data frame.
The code takes time as find_pka(inchikey_string) and find_boiling_point(name) request URL to find the string on the Pubchem page, then extract it using regex. 
The boiling point is a mean of all the values (references) found.
"""
#Finds the pKa using the code of Khoi Van.
def find_pka(inchikey_string):
    text_pka = pka_lookup_pubchem(inchikey_string, "inchikey")
    if text_pka is not None and 'pKa' in text_pka:
        pKa_value = float(text_pka['pKa'])
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


"""
THis code takes the dataframe created with the properties of each compound in the mixture and determines the chromatgraphy from the conditions.
"""
def det_chromato(df):
    if df['Boiling_temp_(°C)'].min() >= 300:
        Chromato_type = 'GC'
        eluent_nature = 'gas'
        proposed_pH = None
    else:
        max_molar_mass = df['Molecular_mass'].max()
        min_pKa = float('inf')  # Initialize min_pKa with a large value
        max_pKa = float('-inf')  # Initialize max_pKa with a small value
        for pKa_entry in df['pKa']:
            if isinstance(pKa_entry, list):
                for pKa_value in pKa_entry:
                    min_pKa = min(pKa_value, min_pKa)
                    max_pKa = max(pKa_value, max_pKa)
            else:
                min_pKa = min(pKa_entry, min_pKa)
                max_pKa = max(pKa_entry, max_pKa)
        
        if max_molar_mass <= 2000:
            max_logP = df['logP'].max()
            min_logP = df['logP'].min()
            if max_logP < 0:
                proposed_pH = max_pKa + 2
                if 3 <= proposed_pH <= 11 and max_pKa + 2 >= proposed_pH:
                    Chromato_type = 'IC'
                    eluent_nature = 'aqueous'
                    proposed_pH = max_pKa + 2
                else:
                    Chromato_type = 'HPLC'
                    eluent_nature = 'organic or hydro-organic'
                    proposed_pH = min_pKa + 2
            else:
                Chromato_type = 'HPLC'
                if -2 <= min_logP <= 0:
                    eluent_nature = 'organic or hydro-organic'
                    if min_logP >= 0:
                        Chromato_type += ' on normal stationary phase'
                    else:
                        Chromato_type += ' on reverse stationary phase using C18 column'
                else:
                    eluent_nature = 'organic or hydro-organic'
                    Chromato_type += ' on normal stationary phase'
                proposed_pH = min_pKa + 2
        else:
            max_logP = df['logP'].max()
            min_logP = df['logP'].min()
            if max_logP < 0:
                Chromato_type = 'HPLC on reverse stationary phase'
                eluent_nature = 'organic or hydro-organic'
                proposed_pH = min_pKa + 2
            else:
                Chromato_type = 'SEC on gel filtration'
                eluent_nature = 'aqueous'
                proposed_pH = min_pKa + 2
    
    return Chromato_type, eluent_nature, proposed_pH




"""
The actual code which is run when calling this file
"""
def main():
    mixture = []
    root = tk.Tk()
    root.title("Determination of Chromatography Type")
    add_molecule()
    
    # Interface widgets
    entry_widget = tk.Entry(root)
    label = tk.Label(root, text="pH value:")
    mixture_entry = ttk.Entry(root)
    mixture_label = ttk.Label(root, text="Names of the molecules in the mixture:")
    add_button = ttk.Button(root, text="Add molecule", command=add_molecule)
    mixture_listbox = tk.Listbox(root)
    calculate_button = ttk.Button(root, text="Determine chromatography", command=lambda: det_chromato(df))
    result1_label = ttk.Label(root, text="")
    result2_label = ttk.Label(root, text="")
    result3_label = ttk.Label(root, text="")
    # Interface layout
    mixture_label.grid(row=0, column=0, padx=5, pady=5)
    mixture_entry.grid(row=0, column=1, padx=5, pady=5)
    add_button.grid(row=1, column=0, columnspan=2, padx=5, pady=5)
    mixture_listbox.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
    calculate_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
    result1_label.grid(row=4, column=0, columnspan=2, padx=5, pady=5)
    result2_label.grid(row=5, column=0, columnspan=2, padx=5, pady=5)
    result3_label.grid(row=6, column=0, columnspan=2, padx=5, pady=5)
    
    root.mainloop()
    
    
    get_df_properties(mixture)
    
    Mixture_chromato_type, eluent_nature, proposed_pH = det_chromato(df)
    print(f"The advisable chromatography type is: {Mixture_chromato_type}")
    print(f"Eluent nature: {eluent_nature}")
    if proposed_pH is not None:
        print(f"Proposed pH for the eluent: {proposed_pH}")


if __name__ == "__main__":
    main()