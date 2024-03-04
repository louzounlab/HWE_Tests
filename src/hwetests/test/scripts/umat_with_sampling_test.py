from hwetests import umat_with_uncertainty

if __name__ == '__main__':
    p_value = umat_with_uncertainty.full_algorithm(file_path='../data/ambiguous_data.csv')
    print(f'p-value: {p_value}')
