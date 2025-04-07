from hwetests import asta
from hwetests.tests import dataloader

if __name__ == '__main__':
    # getting the absolute path to the 'ambiguous_data.csv' file
    ambiguous_data_path = dataloader.get_path(is_ambiguous=True)
    # Perform ASTA
    p_value, statistic, dof = asta.full_algorithm(file_path=ambiguous_data_path,
                                                  cutoff_value=4.0)
    print(f'p-value: {p_value}')
    print(f'statistic: {statistic}')
    print(f'degrees of freedom: {dof}')
