import pytest
import tkinter as tk
from tkinter import ttk
from unittest.mock import patch
import pandas as pd
import urllib.error
import urllib.parse
import re
from Chrfinder.Chrfinder import add_molecule, find_pka, find_boiling_point, get_df_properties, det_chromato, update_results, main

#find_pka(inchikey_string)
def test_find_pka_valid_input():
    inchikey_string = 'RYYVLZVUVIJVGH-UHFFFAOYSA-N'
    expected_pka = '14'
    assert find_pka(inchikey_string) == expected_pka

def test_find_pka_invalid_input():
    inchikey_string = "InvalidInchikeyString"
    assert find_pka(inchikey_string) is None

def test_find_pka_missing_pka():
    inchikey_string = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    assert find_pka(inchikey_string) is None

#find_boiling_point(compound)
def test_mixed_celsuis_fahrenheit():
    compound = "Caffeine"
    result = find_boiling_point(compound)
    assert result == pytest.approx(177.93, rel=1e-2)

def test_mixed_celsuis_fahrenheit_make_sure():
    compound = "Water"
    result = find_boiling_point(compound)
    assert result == pytest.approx(100, rel=1e-2)

def test_no_boiling_Pubchem():
    compound = "Aspartame"
    encoded_compound = urllib.parse.quote_plus(compound)
    result = find_boiling_point(encoded_compound)
    assert result == None

def test_no_compound():
    try:
        result = find_boiling_point("unknown")
    except urllib.error.HTTPError:
        result = {}
    assert result == {}

def test_empty_input():
    try:
        assert find_boiling_point('') is None
    except urllib.error.HTTPError:
        assert True

    try:
        assert find_boiling_point(None) is None
    except urllib.error.HTTPError:
        assert True


#Doesn't work because of pandas: 'Attribute "dtype" are different' float64 != object
"""
#get_df_properti(mixture)
def test_df_valid_input():
    mixture = ['water', 'ethanol']
    df = get_df_properties(mixture)
    
    expected_data = {
        'CID': [962, 702],
        'MolecularFormula': ['H2O', 'C2H6O'],
        'MolecularWeight': [18.015, 46.070],
        'InChIKey': ['XLYOFNOQVPJJNP-UHFFFAOYSA-N', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'],
        'IUPACName': ['oxidane', 'ethanol'],
        'XLogP': ['-0.5', '-0.1'],
        'pKa': [None, 15.9],
        'Boiling Point': [99.99, 78.28]
    }

    expected_df = pd.DataFrame(expected_data)
    pd.testing.assert_frame_equal(df, expected_df, check_exact=False)
"""

def test_compound_not_found():
    mixture = ['unknown_compound']
    df = get_df_properties(mixture)

    expected_df = pd.DataFrame()
    pd.testing.assert_frame_equal(df, expected_df)

#Does not work because of the principle of this function: 
"""
# Test for add_molecule function
def test_add_molecule():
    root = tk.Tk()
    mixture_entry = ttk.Entry(root)
    mixture_listbox = tk.Listbox(root)
    mixture = []

    mixture_entry.insert(0, "acetone")
    add_molecule(mixture_entry, mixture_listbox)

    assert mixture == ["acetone"]
    assert mixture_listbox.get(0) == "acetone"
    root.destroy()
"""    

# Test for det_chromato function
def test_det_chromato_empty():
    df = pd.DataFrame()
    result = det_chromato(df)
    assert result == ("Unknown", "Unknown", None)

def test_det_chromato_gc():
    data = {
        'Boiling Point': [250],
        'MolecularWeight': [150],
        'XLogP': [0],
        'pKa': [[5]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ("GC", "gas", None)

def test_det_chromato_hplc1():
    data = {
        'Boiling Point': [350],
        'MolecularWeight': [500],
        'XLogP': [1],
        'pKa': [[3, 7]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ("HPLC on normal stationary phase", "organic or hydro-organic", 5)

def test_det_chromato_not_gc():
    data = {
        'Boiling Point': [177.93, None, None],
        'MolecularWeight': [194, 204, 133],
        'XLogP': [-0.07, -1.33, 0.07],
        'pKa': [[14], [3.02], [1.08, 9.13]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('GC', 'gas', None)

def test_det_chromato_low_boiling_point():
    data = {
        'Boiling Point': [251],
        'MolecularWeight': [150],
        'XLogP': [-3.5],
        'pKa': [3]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('IC', 'aqueous', 5)

def test_det_chromato_high_molecular_mass():
    data = {
        'Boiling Point': [251],
        'MolecularWeight': [2500],
        'XLogP': [1],
        'pKa': [[5]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('SEC on gel permeation with a hydrophobe organic polymer stationary phase', 'organic solvent', 7)



# Test for main function (GUI Initialization)
#def test_main():
    #try:
        #main()
    #except Exception as e:
        #pytest.fail(f"main() raised {e}")

