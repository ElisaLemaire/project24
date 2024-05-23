<p align="center">
  <img src="assets/Image_Chrfinder.webp" alt="Project Logo" width="650"/>
</p>

# Chrfinder

## <ins>Project overview</ins>

.. image:: http://img.shields.io/pypi/l/PubChemPy.svg?style=flat
    :target: https://github.com/mcs07/PubChemPy/blob/master/LICENSE

Welcome to **Chrfinder**! This project **automates the selection of the most suitable chromatography** technique . By simply providing the **names** of the molecules in the mixture, the code retrieves their physicochemical properties from **PubChem** (web source) and determines the optimal chromatography method based on these properties. It also gives the optimal conditions.

## ✅ <ins>Benefits</ins>

- **🚀 Efficiency**: Automates the property retrieval and decision-making process, saving time and reducing manual effort.
- **🎯 Accuracy**: Utilizes precise physicochemical data to ensure the most suitable chromatography technique is chosen.
- **🌐 Versatility**: Supports a wide range of organic compounds and chromatography methods (PubChem database).


## ⚙ <ins>Installation</ins>

Create a new environment, you may also give the environment a different name. 

```
conda create -n chrfinder python=3.10 
```

```
conda activate chrfinder
pip install .
```

If you need jupyter lab, install it 

```
pip install jupyterlab
```

## 🛠️ <ins>Development installation</ins>

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:Averhv/Chrfinder`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:Averhv/Chrfinder.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
pip install -e ".[test,doc]"
```

#### Run tests and get coverage

```
pip install tox
tox
```

## 📒 <ins>Features</ins>

```python
from chrfinder import main

# Running the whole file ask for molecules through Tkinter and returns the best chromatography
main()
```

#### 🌐 Optional functions

##### find_pka(inchikey)
Finds the pKa value for a compound using its InChIKey.
```python
from chrfinder import find_pka

inchikey = "XEFQLINVKFYRCS-UHFFFAOYSA-N"
find_pka(inchikey)
```

##### find_boiling_point(name)
Finds the boiling point for a compound by name.
```python
from chrfinder import find_boiling_point

compound_name = "Ethanol"
find_boiling_point(compound_name)
```

##### get_df_properties()
Get a DataFrame of properties for a mixture of compounds.
```python
from chrfinder import get_df_properties

mixture = ["Acetone", "Ethanol", "Methanol"]
get_df_properties(mixture, verbose=True)
```

## <ins>How It Works</ins>

1. **Input**: User provides the names of the molecules present in the mixture through a Tkinter interface.

2. **Data Retrieval**: Finds the following key physicochemical properties for each molecule in Pubchem:
     - **Boiling temperature (°C)**
     - **logP (partition coefficient)**
     - **pKa (acid dissociation constants)**
     - **Molecular mass**

3. **Chromatography Type Decision**: Follows logical conditions to determine best chromatography and conditions
   - **Gas Chromatography (GC)**: if the Boiling Point is low (T<sub>eb</sub> <250°C).
   - **Ion Chromatography (IC)**: for small molecules (M<2000g/mol) and a negative maximum LogP negative
     - Selected if the maximum molecular mass is less than or equal to 2000, and the maximum logP is negative, with a proposed pH derived from the pKa values.
   - **High-Performance Liquid Chromatography (HPLC)**: Chosen for different conditions. Stationary phases and eluent natures are suggested.
   - **Size Exclusion Chromatography (SEC)**: For big molecules (M>2000g/mol). From LogP, it suggest gel permeation or gel filtration, with corresponding eluant.

4. **Output**:
   - The code outputs the advisable chromatography type, the nature of the eluent (gas, aqueous, or organic), and the proposed pH for the eluent if applicable through the Tkinter interface.
  
## <ins>Work in progress...</ins>
- Build a data molecules thermostability database;
- taking into account multiple pKa values for polyacids for exemple;
- optimize the research: search only one time te same name;
- find physicalchemical properties as addition functionality;


## 🫱🏽‍🫲🏼 <ins>Contributing</ins>
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

