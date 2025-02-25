import os
import pandas as pd
from hwetests import asta
from hwetests.tests import utils
from matplotlib import pyplot as plt


if __name__ == '__main__':
    data_file_path = '../data/data.freqs'
    directory_path = 'my_directory'

    utils.save_plots_and_data(data_file_path, directory_path)
