from hwetests import umat_with_uncertainty
from hwetests.tests import dataloader

if __name__ == '__main__':
    ambiguous_data_path = dataloader.get_path(
        'ambiguous_data.csv')  # getting the absolute path to the 'ambiguous_data.csv' file
    p_value = umat_with_uncertainty.full_algorithm(file_path='../data/ambiguous_data.csv')
    print(f'p-value: {p_value}')
