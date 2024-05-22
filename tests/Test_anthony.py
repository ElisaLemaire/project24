import sys, os
sys.path.append(os.path.realpath('src'))

import pytest
from unittest.mock import patch
from Chrfinder import find_boiling_point, find_pka
import pandas as pd
from pubchemprops import get_cid_by_name, get_first_layer_props, get_second_layer_props
from pka_lookup import pka_lookup_pubchem


class Test_find_pka:
    @patch('Chrfinder.pka_lookup_pubchem')
    def test_find_pka_valid_input(self, mock_pka_lookup):
        inchikey_string = 'CSCPPACGZOOCGX-UHFFFAOYSA-N'
        expected_pka = '20'
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
            'IUPAC_Name': 'propan-2-one'}
        assert find_pka(inchikey_string) == expected_pka

    @patch('Chrfinder.pka_lookup_pubchem')
    def test_find_pka_invalid_input(self, mock_pka_lookup):
        inchikey_string = "InvalidInchikeyString"
        mock_pka_lookup.return_value = None
        assert find_pka(inchikey_string) == None

    @patch('Chrfinder.pka_lookup_pubchem')
    def test_find_pka_missing_pka(self, mock_pka_lookup):
        inchikey_string = "InchikeyStringMissingPKA"
        mock_pka_lookup.return_value = {
            'source': 'Pubchem',
            'Pubchem_CID': '12345',
            'reference': 'ReferenceInfo',
            'Substance_CASRN': 'CASNumber'
        }
        assert find_pka(inchikey_string) == None



class Test_find_boiling_point:
    def test_single_celsius(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 °C'}]}}]}) == 100
    
    def test_single_fahrenheit(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '212 °F'}]}}]}) == 100
    
    def test_multiple_celsius(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 °C'}]}},
                                                      {'Value': {'StringWithMarkup': [{'String': '50 °C'}]}}]}) == 75
    
    def test_multiple_fahrenheit(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '212 °F'}]}},
                                                      {'Value': {'StringWithMarkup': [{'String': '122 °F'}]}}]}) == 100
    
    def test_mixed_values(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 °C'}]}},
                                                      {'Value': {'StringWithMarkup': [{'String': '212 °F'}]}}]}) == 100
    
    def test_no_data(self):
        assert find_boiling_point({'Some_other_property': [{'Value': {'StringWithMarkup': [{'String': 'random value'}]}}]}) == None
    
    def test_no_value_key(self):
        assert find_boiling_point({'Boiling Point': [{'Key': {'StringWithMarkup': [{'String': '100 °C'}]}}]}) == None
    
    # Invalid Data - No 'StringWithMarkup' Key
    def test_no_string_with_markup_key(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'AnotherKey': [{'String': '100 °C'}]}}]}) == None
    
    # Invalid Data - No Matching Patterns
    def test_no_matching_patterns(self):
        assert find_boiling_point({'Boiling Point': [{'Value': {'StringWithMarkup': [{'String': '100 C'}]}}]}) == None
    
    # Empty Input
    def test_empty_input(self):
        assert find_boiling_point('') == None
        assert find_boiling_point(None) == None



class TestGetDfProperties:
    
    @patch('Chrfinder.get_first_layer_props')
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
        mock_find_pka.side_effect = ['NaN', '15.9']
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
            'pKa': [NaN, 15.9],
            'Boiling Point': [99.99, 78.28]
        }
        expected_df = pd.DataFrame(expected_data)
        pd.testing.assert_frame_equal(df, expected_df)

    
    @patch('Chrfinder.get_first_layer_props')
    def test_compound_not_found(self, mock_get_first_layer_props):
        mock_get_first_layer_props.side_effect = urllib.error.HTTPError(
            url=None, code=404, msg=None, hdrs=None, fp=None
        )
        
        mixture = ['unknown_compound']
        df = get_df_properties(mixture)
        
        expected_df = pd.DataFrame()
        pd.testing.assert_frame_equal(df, expected_df)

    
    @patch('Chrfinder.get_first_layer_props')
    @patch('Chrfinder.find_pka')
    @patch('Chrfinder.find_boiling_point')

    
    def test_no_valid_properties(self, mock_get_first_layer_props):
        mock_get_first_layer_props.return_value = {}
        
        mixture = ['invalid_compound']
        df = get_df_properties(mixture)
        
        expected_df = pd.DataFrame()
        pd.testing.assert_frame_equal(df, expected_df)
        
