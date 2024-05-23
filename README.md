<p align="center">
  <img src="assets/Image_Chrfinder.webp" alt="Project Logo" width="650"/>
</p>

# Chrfinder

## Project overview

Welcome to **Chrfinder**! This project automates the selection of the most suitable chromatography technique for separating and analyzing a given mixture of organic compounds. By simply providing the names of the molecules in the mixture, the code retrieves their physicochemical properties from PubChem web source and determines the optimal chromatography method based on these properties.

## ✅ Benefits

- **🚀 Efficiency**: Automates the property retrieval and decision-making process, saving time and reducing manual effort.
- **🎯 Accuracy**: Utilizes precise physicochemical data to ensure the most suitable chromatography technique is chosen.
- **🌐 Versatility**: Supports a wide range of organic compounds and chromatography methods (PubChem database).

## 🚀 Getting Started

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code).
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## ⚙ Installation

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

## How It Works

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

## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:dsantos03/Chrfinder`.

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

### Run tests and get coverage

```
pip install tox
tox
```

## Work in progress...
- Build a data molecules thermostability database;
- taking into account multiple pKa values for polyacids for exemple;
- optimize the research: search only one time te same name;
- find physicalchemical properties as addition functionality;


## 🫱🏽‍🫲🏼 Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

