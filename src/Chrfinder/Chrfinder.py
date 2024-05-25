import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import pandas as pd
from Chrfinder.pubchemprops import get_cid_by_name, get_first_layer_props, get_second_layer_props
import urllib.error
import urllib.parse
from Chrfinder.pka_lookup import pka_lookup_pubchem
import re
import json

"""
This code takes as input a list of compound written like: acetone, water. The code allows spaces, wrong names and unknown pubchem names.
Then it iterates through each of them to find if they exist on pubchem, and if they do,
then 'CID', 'MolecularFormula', 'MolecularWeight', 'InChIKey', 'IUPACName', 'XLogP', 'pKa',  and 'BoilingPoint' is added into a list and then a data frame.
The code takes time as find_pka(inchikey_string) and find_boiling_point(name) request URL to find the string on the Pubchem page, then extract it using regex. 
The Boiling Point is a mean of all the values (references) found.
"""

mixture=[]

def add_molecule(mixture_entry, mixture_listbox):
    """
    Add a molecule to the mixture list and update the listbox.

    Args:
        mixture_entry (tk.Entry): The entry widget for molecule input.
        mixture_listbox (tk.Listbox): The listbox widget displaying the mixture.

    Returns:
        None
    """
    element = mixture_entry.get()
    mixture.append(element)
    mixture_listbox.insert(tk.END, element.strip())


def add_entry_widget(root):
    """
    Add entry widgets to the root window.

    Args:
        root (tk.Tk): The root window of the Tkinter application.

    Returns:
        None
    """
    entry_widget.grid(row=3, column=1, padx=5, pady=5)
    label.grid(row=3, column=0, padx=5, pady=5)
    
#Finds the pKa using the code of Khoi Van.
def find_pka(inchikey_string, verbose=False):
    """
    Find the pKa value for a compound using its InChIKey.

    Uses the pka_lookup_pubchem() function of the pka_lookup.py file. It converts to float the value (string) of the first pka of pubchem.
    
    Args:
        inchikey_string (str): The InChIKey of the compound.
        verbose (bool, optional): If True, print the pKa value. Defaults to False.

    Returns:
        float or None: The pKa value if found, otherwise None.
    """
    text_pka = pka_lookup_pubchem(inchikey_string, "inchikey")
    if verbose:
        print(text_pka)
    if text_pka is not None and 'pKa' in text_pka:
            pKa_value = text_pka['pKa']
            return pKa_value
    else:
        return None

def find_boiling_point(name, verbose= False):
    """
    Find the boiling point for a compound by name.

    Exctracts the celsius and Farenheit boiling point separetely, using RegExr patterns.
    Converts °F to °C and takes the mean, then returns the float.
    
    Args:
        name (str): The name of the compound.
        verbose (bool, optional): If True, print the boiling point information (used in notebooks). Defaults to False.

    Returns:
        float or None: The average boiling point in Celsius if found, otherwise None.
    """
    text_dict = get_second_layer_props(str(name), ['Boiling Point', 'Vapor Pressure'])
    Boiling_point_values = []
    pattern_celsius = r'([-+]?\d*\.\d+|\d+) °C'
    pattern_F = r'([-+]?\d*\.\d+|\d+) °F'

    if verbose:
        print(text_dict)
    
    if 'Boiling Point' in text_dict:
        for item in text_dict['Boiling Point']:
            if 'Value' in item and 'StringWithMarkup' in item['Value']:
                string_value = item['Value']['StringWithMarkup'][0]['String']
    
                #Search for Celsius values, if found: adds to the list Boiling_point_values
                match_celsius = re.search(pattern_celsius, string_value)
                if match_celsius:
                    celsius = float(match_celsius.group(1))
                    Boiling_point_values.append(celsius)
    
                #Search for Farenheit values, if found: converts farenheit to celsius before adding to the list Boiling_point_values
                match_F = re.search(pattern_F, string_value)
                if match_F:
                    fahrenheit_temp = float(match_F.group(1))
                    celsius_from_F = round(((fahrenheit_temp - 32) * (5/9)), 2)
                    Boiling_point_values.append(celsius_from_F)
                    
        if Boiling_point_values:
            Boiling_temp = round((sum(Boiling_point_values) / len(Boiling_point_values)), 2)
        else:
            Boiling_temp = None
    else:
        Boiling_temp = None
    return Boiling_temp

def get_df_properties(mixture, verbose = False):
    """
    Get a DataFrame of properties for compounds in the mixture.

    Args:
        mixture (list): The list of compound names.
        verbose (bool, optional): If True, print intermediate information (used in notebooks). Defaults to False.

    Returns:
        pd.DataFrame: A DataFrame containing properties of the compounds.
    """
    compound_list = mixture
    compound_properties = []  # Define compound_properties here
    valid_properties = []
    for compound_name in compound_list:
        compound_name_encoded = urllib.parse.quote(compound_name.strip())
        try: 
            first_data = get_first_layer_props(compound_name_encoded, ['MolecularFormula', 'MolecularWeight', 'InChIKey', 'IUPACName', 'XLogP'])
            compound_info = {}
            for prop in ['CID', 'MolecularFormula', 'MolecularWeight', 'InChIKey', 'IUPACName', 'XLogP']:
                if prop == 'MolecularWeight':
                    MolecularWeight_string = first_data.get(prop)
                    if MolecularWeight_string is not None:
                        MolecularWeight_float = float(MolecularWeight_string)
                        compound_info[prop] = MolecularWeight_float
                    else:
                        compound_info[prop] = None
                else:
                    compound_info[prop] = first_data.get(prop)
            if verbose:
                print(compound_info)
            
            #adds pKa if float, else converts to float from string by extracting the float contained (wrongly written on pubchem).
            pka_value = find_pka(first_data['InChIKey'])
            pka_float = None
            if pka_value is not None:
                if isinstance(pka_value, float):
                    pka_float = pka_value
                else:
                    try:
                        pka_float = float(pka_value)
                    except (ValueError, TypeError):
                        match = re.search(r'\d+\.\d+', str(pka_value))
                        if match:
                            pka_float = float(match.group())
            if pka_float is not None:
                compound_info['pKa'] = pka_float
            else:
                # Handle the case where pka_float is None
                compound_info['pKa'] = None
                            
            #adds Boiling point
            compound_info['Boiling Point'] = find_boiling_point(compound_name_encoded)
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
    if verbose:
        print(df)
    return(df)

def det_chromato(df):
    """
    Determine the type of chromatography based on compound properties.

    Conditions: Boiling temperature, molar mass, pka, LogP. To improve: add Thermal stability.

    Args:
        df (pd.DataFrame): DataFrame containing compound properties.

    Returns:
        tuple: Chromatography type, eluent nature, and proposed pH.
    """
    global Type_Label, Eluant_Label, pH_Label
    if df.empty:
        return "Unknown", "Unknown", None
    
    # Filter out NaN values from the boiling points list
    boiling_temps = [temp for temp in df['Boiling Point'] if temp is not None and not pd.isna(temp)] 
    
    # Check if there are valid boiling points and if the maximum is <= 250°C
    if boiling_temps and pd.Series(boiling_temps).max() <= 250:
        Chromato_type = 'GC'
        eluent_nature = 'gas'
        proposed_pH = None
    else:
        molar_masses = [mass for mass in df['MolecularWeight'] if mass is not None and not pd.isna(mass)]
        max_molar_mass = max(molar_masses) if molar_masses else None

        min_pKa = float('inf')
        max_pKa = float('-inf')
        for pKa_entry in df['pKa']:
            if isinstance(pKa_entry, list):
                for pKa_value in pKa_entry:
                    if pKa_value is not None and not pd.isna(pKa_value):
                        min_pKa = min(pKa_value, min_pKa)
                        max_pKa = max(pKa_value, max_pKa)
            else:
                if pKa_entry is not None and not pd.isna(pKa_entry):
                    min_pKa = min(pKa_entry, min_pKa)
                    max_pKa = max(pKa_entry, max_pKa)
        
        if min_pKa == float('inf'):
            min_pKa = None
        if max_pKa == float('-inf'):
            max_pKa = None

        logPs = [XLogP for XLogP in df['XLogP'] if XLogP is not None and not pd.isna(XLogP)]
        max_logP = max(logPs) if logPs else None
        min_logP = min(logPs) if logPs else None

        if max_molar_mass is not None and max_molar_mass <= 2000:
            if max_logP is not None and max_logP < 0:
                proposed_pH = max_pKa + 2 if max_pKa is not None else None
                if proposed_pH is not None and 3 <= proposed_pH <= 11:
                    Chromato_type = 'IC'
                    eluent_nature = 'aqueous'
                else:
                    Chromato_type = 'HPLC'
                    eluent_nature = 'organic or hydro-organic'
                    proposed_pH = min_pKa + 2 if min_pKa is not None else None
            else:
                Chromato_type = 'HPLC'
                if min_logP is not None and -2 <= min_logP <= 0:
                    eluent_nature = 'organic or hydro-organic'
                    if min_logP >= 0:
                        Chromato_type += ' on normal stationary phase'
                    else:
                        Chromato_type += ' on reverse stationary phase using C18 column'
                else:
                    eluent_nature = 'organic or hydro-organic'
                    Chromato_type += ' on normal stationary phase'
                proposed_pH = min_pKa + 2 if min_pKa is not None else None
        else:
            if max_logP is not None and max_logP < 0:
                Chromato_type = 'HPLC on reverse stationary phase'
                eluent_nature = 'organic or hydro-organic'
                proposed_pH = min_pKa + 2 if min_pKa is not None else None
            else:
                if max_logP is not None and max_logP > 0:
                    Chromato_type = 'SEC on gel permeation with a hydrophobe organic polymer stationary phase'
                    eluent_nature = 'organic solvent'
                else:
                    Chromato_type = 'SEC on gel filtration with a polyhydroxylated hydrophile polymer stationary phase'
                    eluent_nature = 'aqueous'
                proposed_pH = min_pKa + 2 if min_pKa is not None else None
    
    return Chromato_type, eluent_nature, proposed_pH

def update_results(root, mixture):
    """
    Update the results displayed in the Tkinter GUI with the determined chromatography type.

    Args:
        root (tk.Tk): The root window of the Tkinter application.
        mixture (list): The list of compounds in the mixture.

    Returns:
        None
    """
    global Type_Label, Eluant_Label, pH_Label
    if not mixture:
        messagebox.showinfo("Error", "Please add molecules to the mixture before determining chromatography.")
        return
    
    df = get_df_properties(mixture)
    Chromato_type, eluent_nature, proposed_pH = det_chromato(df)
    
    Type_Label.config(text=f"The advisable chromatography type is: {Chromato_type}")
    Eluant_Label.config(text=f"Eluent nature: {eluent_nature}")
    if proposed_pH is not None:
        pH_Label.config(text=f"Proposed pH for the eluent: {proposed_pH}")

def main():
    """
    The main function to initialize and run the Tkinter GUI application.

    Returns:
        None
    """
    global entry_widget, label, mixture_listbox, Type_Label, Eluant_Label, pH_Label
    root = tk.Tk()
    root.title("Determination of Chromatography Type")
    """
    get_df_properties(mixture_test)
    Mixture_chromato_type, eluent_nature, proposed_pH = det_chromato(df)
    """
    entry_widget = tk.Entry(root)
    label = tk.Label(root, text="pH value:")
    mixture_entry = ttk.Entry(root)
    mixture_label = ttk.Label(root, text="Names of the molecules in the mixture:")
    add_button = ttk.Button(root, text="Add molecule", command=lambda: add_molecule(mixture_entry, mixture_listbox))
    mixture_listbox = tk.Listbox(root)
    calculate_button = ttk.Button(root, text="Determine chromatography", command=lambda: update_results(root, mixture))
    Type_Label = ttk.Label(root, text="")
    Eluant_Label = ttk.Label(root, text="")
    pH_Label = ttk.Label(root, text="")
    
    mixture_label.grid(row=0, column=0, padx=5, pady=5)
    mixture_entry.grid(row=0, column=1, padx=5, pady=5)
    add_button.grid(row=1, column=0, columnspan=2, padx=5, pady=5)
    mixture_listbox.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
    calculate_button.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
    Type_Label.grid(row=4, column=0, columnspan=2, padx=5, pady=5)
    Eluant_Label.grid(row=5, column=0, columnspan=2, padx=5, pady=5)
    pH_Label.grid(row=6, column=0, columnspan=2, padx=5, pady=5)
    
    root.mainloop()

if __name__ == "__main__":
    main()
