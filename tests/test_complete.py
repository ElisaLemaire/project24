import pytest
import tkinter as tk
from tkinter import ttk
from unittest.mock import patch
import pandas as pd
import urllib.error
import urllib.parse
import re
import sys, os
sys.path.append(os.path.realpath('src'))
from Chrfinder import add_molecule, find_pka, find_boiling_point, get_df_properties, det_chromato, update_results, main

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

# Test for find_pka function
@patch('main_script.pka_lookup_pubchem')
def test_find_pka_valid(mock_pka_lookup_pubchem):
    mock_pka_lookup_pubchem.return_value = {'pKa': '7.4'}
    assert find_pka("CSCPPACGZOOCGX-UHFFFAOYSA-N") == '7.4'

@patch('main_script.pka_lookup_pubchem')
def test_find_pka_invalid(mock_pka_lookup_pubchem):
    mock_pka_lookup_pubchem.return_value = None
    assert find_pka("InvalidInChIKey") is None

# Test for find_boiling_point function
@patch('main_script.get_second_layer_props')
def test_find_boiling_point(mock_get_second_layer_props):
    mock_get_second_layer_props.return_value = {
        'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 °C'}]}}]
    }
    assert find_boiling_point("water") == 100

@patch('main_script.get_second_layer_props')
def test_find_boiling_point_multiple(mock_get_second_layer_props):
    mock_get_second_layer_props.return_value = {
        'Boiling Point': [
            {'Value': {'StringWithMarkup': [{'String': '100 °C'}]}},
            {'Value': {'StringWithMarkup': [{'String': '212 °F'}]}}
        ]
    }
    assert find_boiling_point("water") == 100

@patch('main_script.get_second_layer_props')
def test_find_boiling_point_none(mock_get_second_layer_props):
    mock_get_second_layer_props.return_value = {}
    assert find_boiling_point("unknown") is None

# Test for get_df_properties function
@patch('main_script.get_first_layer_props')
@patch('main_script.find_pka')
@patch('main_script.find_boiling_point')
def test_get_df_properties(mock_find_boiling_point, mock_find_pka, mock_get_first_layer_props):
    mock_get_first_layer_props.return_value = {
        'CID': '962',
        'MolecularFormula': 'H2O',
        'MolecularWeight': '18.015',
        'InChIKey': 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
        'IUPACName': 'oxidane',
        'XLogP': '-0.5'
    }
    mock_find_pka.return_value = '15.9'
    mock_find_boiling_point.return_value = 100

    mixture = ['water']
    df = get_df_properties(mixture)
    
    expected_data = {
        'CID': ['962'],
        'MolecularFormula': ['H2O'],
        'MolecularWeight': [18.015],
        'InChIKey': ['XLYOFNOQVPJJNP-UHFFFAOYSA-N'],
        'IUPACName': ['oxidane'],
        'XLogP': ['-0.5'],
        'pKa': [15.9],
        'Boiling Point': [100]
    }
    expected_df = pd.DataFrame(expected_data)
    pd.testing.assert_frame_equal(df, expected_df)

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

def test_det_chromato_hplc():
    data = {
        'Boiling Point': [350],
        'MolecularWeight': [500],
        'XLogP': [1],
        'pKa': [[3, 7]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ("HPLC on reverse stationary phase using C18 column", "organic or hydro-organic", 5)

# Test for update_results function
@patch('main_script.get_df_properties')
@patch('main_script.det_chromato')
def test_update_results(mock_det_chromato, mock_get_df_properties):
    root = tk.Tk()
    mock_get_df_properties.return_value = pd.DataFrame()
    mock_det_chromato.return_value = ("HPLC", "organic or hydro-organic", 5)
    
    mixture = ["water"]
    update_results(root, mixture)
    
    assert "HPLC" in root.winfo_children()[4].cget("text")
    assert "organic or hydro-organic" in root.winfo_children()[5].cget("text")
    assert "5" in root.winfo_children()[6].cget("text")
    root.destroy()

# Test for main function (GUI Initialization)
def test_main():
    try:
        main()
    except Exception as e:
        pytest.fail(f"main() raised {e}")
