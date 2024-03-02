from hwetests import asta

if __name__ == '__main__':
    p_value, statistic, dof = asta.full_algorithm(file_path='../data/ambiguous_data.csv',
                                                  cutoff_value=4.0)
    print(f'p-value: {p_value}')
    print(f'statistic: {statistic}')
    print(f'degrees of freedom: {dof}')
