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

-  [Usage](#usage)
- [Tests](#tests)

-  [Credits](#credits)

-  [License](#license)

## Installation
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install hwetests.
```bash
pip install hwetests
```

## How to use our package?
Our package is very simple to use. 
in order to run any statistical test, you only need to use one function to get your results.
One can find the statistical tests in

Here is the tour for each statistical test:
### ASTA


## Demos

Here we show how to use our package with simulated data.
You can find the scripts and the simulated data in:
```bash
├───src
│   ├───test
│   │   ├───data
│   │   │   └───unambiguous_data.csv # for ASTA and UMAT with sampling
│   │   │   └───ambiguous_data.csv # for UMAT
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

## Tests

Go the extra mile and write tests for your application. Then provide examples on how to run them here.

## Badges
![Static Badge](https://img.shields.io/badge/python-100%25-blue)


## License
asasda