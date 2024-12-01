import os
import pandas as pd
from HWE_Tests.hwetests import asta
from HWE_Tests.hwetests.tests import utils
from matplotlib import pyplot as plt


if __name__ == '__main__':
    data_file_path = '../data/data.freqs'
    directory_path = 'my_directory'

    utils.save_plots_and_data(data_file_path, directory_path)
