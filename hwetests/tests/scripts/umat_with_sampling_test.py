from HWE_Tests.hwetests import umat_with_uncertainty
from HWE_Tests.hwetests.tests import dataloader

if __name__ == '__main__':
    # getting the absolute path to the 'ambiguous_data.csv' file
    ambiguous_data_path = dataloader.get_path(is_ambiguous=True)
    # Perform UMAT with sampling
    p_value = umat_with_uncertainty.full_algorithm(file_path='../data/ambiguous_data.csv')
    print(f'p-value: {p_value}')
