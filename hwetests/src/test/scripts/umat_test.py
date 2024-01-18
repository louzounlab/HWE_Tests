from hwetests import umat
import numpy as np

if __name__ == '__main__':
    # import data from csv file as a numpy array
    data = np.genfromtxt('../data/unambiguous_data.csv', delimiter=',')
    # run the test
    p_value = umat.full_algorithm(data)
    print(f'p-value: {p_value}')
