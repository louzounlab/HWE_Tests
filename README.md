<p align="center">
    <br>
    <img src="https://raw.githubusercontent.com/huggingface/diffusers/main/docs/source/en/imgs/diffusers_library.jpg" width="400"/>
    <br>
<p>
<p align="center">
<p align="center">
    <a href="https://img.shields.io/badge/python-100%25-blue">
        <img alt="python" src="https://img.shields.io/badge/python-100%25-blue">
    </a>
    <a href="https://img.shields.io/badge/license-MIT-blue">
        <img alt="license" src="https://img.shields.io/badge/license-MIT-blue">
    </a>


# hwetests
![Static Badge](https://img.shields.io/badge/python-100%25-blue) ![Static Badge](https://img.shields.io/badge/license-MIT-blue
)

The Hardy-Weinberg Equilibrium (HWE) assumption is essential to many population genetics models, which assumes that allele pairs from a given population are random. An HWE test needs to test whether the pairing are random or not.

Our python package contains three statistical tests for HWE testing:
- **ASTA**
- **UMAT**
- **UMAT with sampling**

Both  **ASTA** and **UMAT with sampling** assume ambiguity in the observations while **UMAT** does not.

## Table of Contents

-  [Installation](#installation)

-  [Quick tour](#quick_tour)
- [Examples](#examples)

-  [Credits](#credits)

-  [License](#license)

## Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install hwetests.
```bash
pip install hwetests
```

## Quick tour
To immediately use our package, you only need to run a single function.

### ASTA
```python
def full_algorithm(file_path,  
                   is_first_row_contains_columns_names=False,  
                   cutoff_value=0.0,  
                   should_save_csv=False,  
                   should_save_plot=False,  
                   title=''):
```
Suppose

## Examples

Here we show how to use our package with simulated data.
You can find the scripts and the simulated data in:
```bash
├───src
│   ├───test
│   │   ├───data
│   │   │   └───unambiguous_data.csv # for ASTA and UMAT with sampling (contains 50k population, 20 alleles, 0.2 uncertainty, in HWE)
│   │   │   └───ambiguous_data.csv # for UMAT (contains 100k population, in HWE)
│   │   ├───scripts
│   │   │   └───asta_test.py
│   │   │   └───umat_test.py
│   │   │   └───umat_with_sampling_test.py
```

### ASTA test
```python
from hwetests import asta  
  
if __name__ == '__main__':  
    p_value, statistic, dof = asta.full_algorithm(file_path='../data/ambiguous_data.csv',  
                                                  cutoff_value=4.0)  
    print(f'p-value: {p_value}')  
    print(f'statistic: {statistic}')  
    print(f'degrees of freedom: {dof}')
```

### UMAT test
```python
from hwetests import umat  
import numpy as np  
  
if __name__ == '__main__':  
    # import data from csv file as a numpy array  
    data = np.genfromtxt('../data/unambiguous_data.csv', delimiter=',')  
    # run the test  
    p_value = umat.full_algorithm(data)  
    print(f'p-value: {p_value}')
```

### UMAT with sampling test
```python
from hwetests import umat_with_uncertainty  
  
if __name__ == '__main__':  
    p_value = umat_with_uncertainty.full_algorithm(file_path='../data/ambiguous_data.csv')
```


## License
asasda
