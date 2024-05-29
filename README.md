<p align="center">
  <img src="assets/Image_Chrfinder.webp" alt="Project Logo" width="650"/>
</p>

# Chrfinder

## <ins>Project overview</ins>

<h1 align="center">
    
[![PyPI version](https://img.shields.io/pypi/v/Chrfinder?style=plastic&color=blue)](https://pypi.python.org/pypi/Chrfinder) 
[![License](https://img.shields.io/github/license/Averhv/Chrfinder?style=plastic&color=Orange)](https://github.com/Averhv/Chrfinder/blob/master/LICENSE) 
<img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/Chrfinder?style=plastic" /> 
[![Coverage Badge](https://github.com/Averhv/Chrfinder/blob/main/assets/coverage-badge.svg?short_path=eb9c651)](https://github.com/Averhv/Chrfinder)


</h1>

Welcome to **Chrfinder**! This project **automates the selection of the most suitable chromatography** technique . By simply providing the **names** of the molecules in the mixture, the code retrieves their physicochemical properties from **PubChem** (web source) and determines the optimal chromatography method based on these properties. It also gives the optimal conditions.

## ✅ <ins>Benefits</ins>

- **🚀 Efficiency**: Automates the property retrieval and decision-making process, saving time and reducing manual effort.
- **🎯 Accuracy**: Utilizes precise physicochemical data to ensure the most suitable chromatography technique is chosen.
- **🌐 Versatility**: Supports a wide range of organic compounds and chromatography methods (PubChem database).


## ⚙ <ins>Installation</ins>

To get started with Chrfinder, you can follow these steps:
Create a new environment, you may also give the environment a different name. 

```bash
conda create -n Chrfinder python=3.10 
```
Then install it through pypi (easy, recommended):
```bash
conda activate Chrfinder
```
```bash
pip install Chrfinder
```
You also have the choice to install it without pypi:

```bash
git clone https://github.com/Averhv/Chrfinder.git
cd Chrfinder
pip install .
```

## 🛠️ <ins>Installation for Development</ins>

If you want to contribute to Chrfinder or run the tests and get coverage, you can install the package in editable mode along with the necessary dependencies for testing and documentation. The following script allows the changes to be reflected immediately:
```bash
git clone https://github.com/Averhv/Chrfinder.git
cd Chrfinder
pip install -e ".[test,doc]"
```

Then you need to run the tests as follow in your terminal:
```bash
pip install tox
tox
```
Test result: 15 passed in ~20s
## 📒 <ins>Features</ins>

The following should be written in the terminal in python mode:
```python
from Chrfinder import main
```
Running the main file asks for molecules through Tkinter and returns the best chromatography.
```python
main()
```

#### 🌐 Optional functions

- find_pka(inchikey)
Finds the pKa value for a compound using its InChIKey.
```python
from Chrfinder import find_pka

inchikey = "XEFQLINVKFYRCS-UHFFFAOYSA-N"
find_pka(inchikey)
```

- find_boiling_point(name)
Finds the boiling point for a compound by name.
```python
from Chrfinder import find_boiling_point

compound_name = "Ethanol"
find_boiling_point(compound_name)
```

- get_df_properties()
Get a DataFrame of properties for a mixture of compounds.
```python
from Chrfinder import get_df_properties

mixture = ["Acetone", "Ethanol", "Methanol"]
get_df_properties(mixture, verbose=True)
```
  
## <ins>Work in progress...</ins>
- **Efficiency: build a database**;
- **Research thermostability** in lab to improve precision
- Taking into account **multiple pKa values** for polyacids for exemple;
- **Optimize the research**: search only one time the same name;


## 🫱🏽‍🫲🏼 <ins>Contributing</ins>
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

