# Chrfinder

<br>

## Project overview

This project automates the selection of the most suitable chromatography technique for analyzing a given mixture of organic compounds. By simply providing the names of the molecules in the mixture, the code retrieves their physicochemical properties from various web sources and determines the optimal chromatography method based on these properties.

## Benefits

- **Efficiency**: Automates the property retrieval and decision-making process, saving time and reducing manual effort.
- **Accuracy**: Utilizes precise physicochemical data to ensure the most suitable chromatography technique is chosen.
- **Versatility**: Supports a wide range of organic compounds and chromatography methods.


## 🔥 Getting Started

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code).
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n chrfinder python=3.10 
```

```
conda activate chrfinder
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(chrfinder) $ pip install jupyterlab
```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:dsantos03/Chrfinder`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:dsantos03/Chrfinder.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(chrfinder) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```

## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.

