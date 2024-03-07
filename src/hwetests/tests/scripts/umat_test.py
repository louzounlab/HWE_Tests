from hwetests import umat
from hwetests.tests import dataloader
import numpy as np

if __name__ == '__main__':
    unambiguous_data_path = dataloader.get_path('unambiguous_data.csv')  # getting the absolute path to the 'unambiguous_data.csv' file
    # import data from csv file as a numpy array
    data = np.genfromtxt(unambiguous_data_path, delimiter=',')
    # run the tests
    p_value = umat.full_algorithm(data)
    print(f'p-value: {p_value}')
