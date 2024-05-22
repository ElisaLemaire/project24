import pandas as pd
import pytest
from src.chrfinder.chromatography import det_chromato

@pytest.mark.parametrize("description, input_data, expected_output", [
    (
        'Basic Functionality Test',
        {
            'Mixture': ['Caf', 'Ace', 'Asp'],
            'Boiling_temp_(°C)': [178, 332.7, 246],
            'logP': [-0.07, -1.33, 0.07],
            'pKa': [[14], [3.02], [1.08, 9.13]],
            'Molecular_mass': [194, 204, 133],
        },
        ('HPLC', 'organic or hydro-organic', 3.08)
    ),
    (
        'High Boiling Point Test',
        {
            'Mixture': ['X'],
            'Boiling_temp_(°C)': [310],
            'logP': [0],
            'pKa': [[5]],
            'Molecular_mass': [150],
        },
        ('GC', 'gas', None)
    ),
    (
        'Low Molecular Mass Test',
        {
            'Mixture': ['Y'],
            'Boiling_temp_(°C)': [100],
            'logP': [-0.5],
            'pKa': [[3]],
            'Molecular_mass': [150],
        },
        ('IC', 'aqueous', 5)
    ),
    (
        'Mixed pKa Values Test',
        {
            'Mixture': ['Z'],
            'Boiling_temp_(°C)': [200],
            'logP': [0],
            'pKa': [[2, 6]],
            'Molecular_mass': [180],
        },
        ('HPLC', 'organic or hydro-organic', 4)
    ),
    (
        'Edge Cases for pKa and logP',
        {
            'Mixture': ['W'],
            'Boiling_temp_(°C)': [250],
            'logP': [-2],
            'pKa': [[3, 11]],
            'Molecular_mass': [160],
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
            'Boiling_temp_(°C)': [None],
            'logP': [None],
            'pKa': [[None]],
            'Molecular_mass': [None],
        },
        (None, None, None)
    ),
    (
        'High Molecular Mass Test',
        {
            'Mixture': ['U'],
            'Boiling_temp_(°C)': [150],
            'logP': [1],
            'pKa': [[5]],
            'Molecular_mass': [2500],
        },
        ('SEC on gel permeation with a hydrophobe organic polymer stationary phase', 'organic solvent', 7)
    ),
])
def test_det_chromato(description, input_data, expected_output):
    df = pd.DataFrame(input_data)
    result = det_chromato(df)
    assert result == expected_output, f"Test failed: {description}.\nExpected {expected_output}, got {result}"
