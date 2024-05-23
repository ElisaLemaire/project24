import pytest
import tkinter as tk
from tkinter import ttk
from unittest.mock import patch
import pandas as pd
import urllib.error
import urllib.parse
import re
from Chrfinder.Chrfinder import add_molecule, find_pka, find_boiling_point, get_df_properties, det_chromato, update_results, main

class TestFindPka:
    @patch('Chrfinder.pka_lookup.pka_lookup_pubchem')
    def test_find_pka_valid_input(self, mock_pka_lookup):
        inchikey_string = 'CSCPPACGZOOCGX-UHFFFAOYSA-N'
        expected_pka = '20.0'
        mock_pka_lookup.return_value = {
            'source': 'Pubchem',
            'Pubchem_CID': '180',
            'pKa': '20',
            'reference': 'Serjeant EP, Dempsey B; Ionisation constants of organic acids in aqueous solution. IUPAC Chem Data Ser No.23. NY,NY: Pergamon pp.989 (1979)',
            'Substance_CASRN': '67-64-1',
            'Canonical_SMILES': 'CC(=O)C',
            'Isomeric_SMILES': 'CC(=O)C',
            'InChI': 'InChI=1S/C3H6O/c1-3(2)4/h1-2H3',
            'InChIKey': 'CSCPPACGZOOCGX-UHFFFAOYSA-N',
            'IUPAC_Name': 'propan-2-one'
        }
        assert find_pka(inchikey_string) == expected_pka

    @patch('Chrfinder.pka_lookup.pka_lookup_pubchem')
    def test_find_pka_invalid_input(self, mock_pka_lookup):
        inchikey_string = "InvalidInchikeyString"
        mock_pka_lookup.return_value = None
        assert find_pka(inchikey_string) is None

    @patch('Chrfinder.pka_lookup.pka_lookup_pubchem')
    def test_find_pka_missing_pka(self, mock_pka_lookup):
        inchikey_string = "InchikeyStringMissingPKA"
        mock_pka_lookup.return_value = {
            'source': 'Pubchem',
            'Pubchem_CID': '12345',
            'reference': 'ReferenceInfo',
            'Substance_CASRN': 'CASNumber'
        }
        assert find_pka(inchikey_string) is None


class TestFindBoilingPoint:
    @patch('Chrfinder.get_second_layer_props')
    def test_single_celsius(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {
            'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 °C'}]}}]
        }
        assert find_boiling_point("water") == 100

    @patch('Chrfinder.get_second_layer_props')
    def test_single_fahrenheit(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {
            'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '212 °F'}]}}]
        }
        assert find_boiling_point("water") == 100

    @patch('Chrfinder.get_second_layer_props')
    def test_multiple_celsius(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {
            'Boiling Point': [
                {'Value': {'StringWithMarkup': [{'String': '100 °C'}]}},
                {'Value': {'StringWithMarkup': [{'String': '50 °C'}]}}
            ]
        }
        assert find_boiling_point("water") == 75

    @patch('Chrfinder.get_second_layer_props')
    def test_multiple_fahrenheit(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {
            'Boiling Point': [
                {'Value': {'StringWithMarkup': [{'String': '212 °F'}]}},
                {'Value': {'StringWithMarkup': [{'String': '122 °F'}]}}
            ]
        }
        assert find_boiling_point("water") == 100

    @patch('Chrfinder.get_second_layer_props')
    def test_mixed_values(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {
            'Boiling Point': [
                {'Value': {'StringWithMarkup': [{'String': '100 °C'}]}},
                {'Value': {'StringWithMarkup': [{'String': '212 °F'}]}}
            ]
        }
        assert find_boiling_point("water") == 100

    @patch('Chrfinder.get_second_layer_props')
    def test_no_data(self, mock_get_second_layer_props):
        mock_get_second_layer_props.return_value = {}
        assert find_boiling_point("unknown") is None

    def test_no_value_key(self):
        assert find_boiling_point({'Boiling Point': [{'Key': {'StringWithMarkup': [{'String': '100 °C'}]}}]}) is None

    def test_no_string_with_markup_key(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'AnotherKey': [{'String': '100 °C'}]}}]}) is None

    def test_no_matching_patterns(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 C'}]}}]}) is None

    def test_empty_input(self):
        assert find_boiling_point('') is None
        assert find_boiling_point(None) is None


class TestGetDfProperties:
    @patch('Chrfinder.pubchemprops.get_first_layer_props')
    @patch('Chrfinder.find_pka')
    @patch('Chrfinder.find_boiling_point')
    def test_valid_input(self, mock_find_boiling_point, mock_find_pka, mock_get_first_layer_props):
        mock_get_first_layer_props.side_effect = [
            {
                'CID': '962',
                'MolecularFormula': 'H2O',
                'MolecularWeight': '18.015',
                'InChIKey': 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
                'IUPACName': 'oxidane',
                'XLogP': '-0.5'
            },
            {
                'CID': '702',
                'MolecularFormula': 'C2H6O',
                'MolecularWeight': '46.070',
                'InChIKey': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
                'IUPACName': 'ethanol',
                'XLogP': '-0.1'
            }
        ]
        mock_find_pka.side_effect = [None, '15.9']
        mock_find_boiling_point.side_effect = [99.99, 78.28]

        mixture = ['water', 'ethanol']
        df = get_df_properties(mixture)

        expected_data = {
            'CID': ['962', '702'],
            'MolecularFormula': ['H2O', 'C2H6O'],
            'MolecularWeight': [18.015, 46.070],
            'InChIKey': ['XLYOFNOQVPJJNP-UHFFFAOYSA-N', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'],
            'IUPACName': ['oxidane', 'ethanol'],
            'XLogP': ['-0.5', '-0.1'],
            'pKa': [None, 15.9],
            'Boiling Point': [99.99, 78.28]
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(df, expected_df)

    @patch('Chrfinder.pubchemprops.get_first_layer_props')
    def test_compound_not_found(self, mock_get_first_layer_props):
        mock_get_first_layer_props.side_effect = urllib.error.HTTPError(
            url=None, code=404, msg=None, hdrs=None, fp=None
        )

        mixture = ['unknown_compound']
        df = get_df_properties(mixture)

        expected_df = pd.DataFrame()
        pd.testing.assert_frame_equal(df, expected_df)

    @patch('Chrfinder.pubchemprops.get_first_layer_props')
    def test_no_valid_properties(self, mock_get_first_layer_props):
        mock_get_first_layer_props.return_value = {}

        mixture = ['invalid_compound']
        df = get_df_properties(mixture)

        expected_df = pd.DataFrame()
        pd.testing.assert_frame_equal(df, expected_df)


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
    assert result == ("HPLC on reverse stationary phase using C18 column", "organic or hydro-organic", 5)

def test_det_chromato_hplc2():
    data = {
        'Boiling Point': [178, 332.7, 246],
        'MolecularWeight': [194, 204, 133],
        'XLogP': [-0.07, -1.33, 0.07],
        'pKa': [[14], [3.02], [1.08, 9.13]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('HPLC', 'organic or hydro-organic', 3.08)

def test_det_chromato_low_boiling_point():
    data = {
        'Boiling Point': [100],
        'MolecularWeight': [150],
        'XLogP': [-0.5],
        'pKa': [3]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('IC', 'aqueous', 5)

def test_det_chromato_high_molecular_mass():
    data = {
        'Boiling Point': [150],
        'MolecularWeight': [2500],
        'XLogP': [1],
        'pKa': [[5]]
    }
    df = pd.DataFrame(data)
    result = det_chromato(df)
    assert result == ('SEC on gel permeation with a hydrophobe organic polymer stationary phase', 'organic solvent', 7)

# Test for update_results function
@patch('Chrfinder.get_df_properties')
@patch('Chrfinder.det_chromato')
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
#def test_main():
    #try:
        #main()
    #except Exception as e:
        #pytest.fail(f"main() raised {e}")

