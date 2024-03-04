from hwetests import asta
from hwetests.test import dataloader

if __name__ == '__main__':
    ambiguous_data_path = dataloader.get_path('ambiguous_data.csv')
    p_value, statistic, dof = asta.full_algorithm(file_path=ambiguous_data_path,
                                                  cutoff_value=4.0)
    print(f'p-value: {p_value}')
    print(f'statistic: {statistic}')
    print(f'degrees of freedom: {dof}')