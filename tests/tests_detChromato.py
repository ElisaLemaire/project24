import pandas as pd
import pytest
import sys, os
sys.path.append(os.path.realpath('src'))
from unittest.mock import patch

class TestDetChromato:
    @pytest.mark.parametrize("description, input_data, expected_output", [
        (
            'Basic Functionality Test',
            {
                'Mixture': ['Caf', 'Ace', 'Asp'],
                'Boiling Point': [178, 332.7, 246],
                'XlogP': [-0.07, -1.33, 0.07],
                'pKa': [[14], [3.02], [1.08, 9.13]],
                'MolecularWeight': [194, 204, 133],
            },
            ('HPLC', 'organic or hydro-organic', 3.08)
        ),
        (
            'High Boiling Point Test',
            {
                'Mixture': ['X'],
                'Boiling Point': [310],
                'XlogP': [0],
                'pKa': [[5]],
                'MolecularWeight': [150],
            },
            ('GC', 'gas', None)
        ),
        (
            'Low Molecular Mass Test',
            {
                'Mixture': ['Y'],
                'Boiling Point': [100],
                'XlogP': [-0.5],
                'pKa': [[3]],
                'MolecularWeight': [150],
            },
            ('IC', 'aqueous', 5)
        ),
        (
            'Mixed pKa Values Test',
            {
                'Mixture': ['Z'],
                'Boiling Point': [200],
                'XlogP': [0],
                'pKa': [[2, 6]],
                'MolecularWeight': [180],
            },
            ('HPLC', 'organic or hydro-organic', 4)
        ),
        (
            'Edge Cases for pKa and logP',
            {
                'Mixture': ['W'],
                'Boiling Point': [250],
                'XlogP': [-2],
                'pKa': [[3, 11]],
                'MolecularWeight': [160],
            },
            ('HPLC on reverse stationary phase using C18 column', 'organic or hydro-organic', 5)
        ),
        (
            'Empty DataFrame Test',
            {},
            (None, None, None)
        ),
        (
            'None Values in DataFrame',
            {
                'Mixture': ['V'],
                'Boiling Point': [None],
                'XlogP': [None],
                'pKa': [[None]],
                'MolecularWeight': [None],
            },
            (None, None, None)
        ),
        (
            'High Molecular Mass Test',
            {
                'Mixture': ['U'],
                'Boiling Point': [150],
                'XlogP': [1],
                'pKa': [[5]],
                'MolecularWeight': [2500],
            },
            ('SEC on gel permeation with a hydrophobe organic polymer stationary phase', 'organic solvent', 7)
        ),
    ])
    def test_det_chromato(self, description, input_data, expected_output):
        df = pd.DataFrame(input_data)
        result = det_chromato(df)
        assert result == expected_output, f"Test failed: {description}.\nExpected {expected_output}, got {result}"

        
