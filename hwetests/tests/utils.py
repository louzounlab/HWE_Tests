import os
import pandas as pd
from HWE_Tests.hwetests import asta
from matplotlib import pyplot as plt

################## CONSTANTS ##################
# levels of the alleles
LEVELS_LIST = [
    'A',
    'B',
    'C',
    'DQB1',
    'DRB1'
]

# dictionary of levels -> relative positions in the real data file
LEVELS_DICT = {'A': 0,
               'B': 1,
               'C': 2,
               'DQB1': 3,
               'DRB1': 4}

RACES_1_LIST = ['ALL']


################## END OF CONSTANTS ##################


def _create_directories(directory_path=''):
    # path to this directory
    # current_path = os.getcwd()

    # -> make directory: data -> real_data -> levels
    # levels_path = os.path.join(current_path, "data", "sapir_data", directory_name, "levels")
    levels_path = os.path.join(directory_path, "levels")
    if not os.path.exists(levels_path):
        os.makedirs(levels_path)

    # get levels, races
    levels_list = LEVELS_LIST
    races_list = RACES_1_LIST

    # create data -> real_data_ races -> {race}
    races_path = os.path.join(directory_path, "races")
    if not os.path.exists(races_path):
        os.makedirs(races_path)
    for race in races_list:
        current_race_path = os.path.join(races_path, race)
        if not os.path.exists(current_race_path):
            os.makedirs(current_race_path)

    for level in levels_list:
        # for each level create directory: levels -> level
        current_level_path = os.path.join(levels_path, level)
        if not os.path.exists(current_level_path):
            os.makedirs(current_level_path)

        # -> races
        races_path = os.path.join(current_level_path, "races")
        # for each level and race create directory: levels -> level -> races -> race
        for race in races_list:
            current_race_path = os.path.join(races_path, race)
            if not os.path.exists(current_race_path):
                os.makedirs(current_race_path)

    results_path = os.path.join(directory_path, "results")
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    plots_path = os.path.join(directory_path, "plots")
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)


def _create_sums_files(data_file_path, directory_path=''):  # freqs file
    races_list = RACES_1_LIST

    # data_file_path = f'../data/sapir_data/{directory_name}.freqs'

    race_to_id_sum = {}

    current_race = ''
    current_id = ''
    current_sum = 0.0

    # races_df = pd.read_csv('data/id_race_don_us.csv',
    #                        dtype={'id': str, 'rcae': str})

    with open(data_file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 500 == 0:
                print(index)
            lst = line.strip('\n').split(',')

            id_row = lst[0]
            observed_row = float(lst[2])
            person_index_row = int(lst[3])

            # if we got to a new id
            if person_index_row == 0 and current_id:
                # save the last person
                if current_race not in race_to_id_sum:
                    race_to_id_sum[current_race] = []
                # race_to_id_amount[current_race][current_id] = current_amount
                race_to_id_sum[current_race].append({'id': current_id, 'sum': current_sum})
                current_sum = 0.0

            current_race = _get_race(id_row)
            current_id = id_row
            current_sum += observed_row

    # save also the last person
    if current_race not in race_to_id_sum:
        race_to_id_sum[current_race] = []
    race_to_id_sum[current_race].append({'id': current_id, 'sum': current_sum})

    # save the files
    for race in races_list:
        list_of_dicts = race_to_id_sum[race]
        df = pd.DataFrame.from_dict(list_of_dicts)
        path_save = os.path.join(directory_path, 'races', race, 'id_sum')
        df.to_csv(path_save, index=False)


# every person appears in subsequent rows with an index indicating how many times he appears.
# for every level and race we create a probabilities file that contains ID,ALLELE_1,ALLELE_2,O_K(I,J)
def _create_id_allele1_allele2_probability_files_for_each_level(level_,
                                                                data_file_path,
                                                                directory_path=''):
    # dictionary (race -> dict of probabilities [ID,ALLELE_1,ALLELE_2,O_K(I,J)] )
    # race_to_list_of_dicts = {}

    # races_df = pd.read_csv('data/id_race_don_us.csv',
    #                        dtype={'id': str, 'rcae': str})
    races_list = RACES_1_LIST

    # data_file_path = f'../data/sapir_data/{directory_name}.freqs'
    # directory_path = f'../data/sapir_data/{directory_name}'

    # race -> dictionary (id -> sum)
    print(f'Creating id amount dictionary:')
    race_to_id_to_sum = {}
    for race in races_list:
        # dtypes = {'id': str, 'sum': float}
        # df = pd.read_csv(f'{real_data_path}/races/{race}/id_sum',
        #                  dtype=dtypes)
        my_dict = {}
        with open(os.path.join(directory_path, 'races', race, 'id_sum'), encoding="utf8") as infile:
            for index, line in enumerate(infile):
                if index == 0:
                    continue
                lst = line.strip('\n').split(',')
                id_row = lst[0]
                sum_row = lst[1]
                my_dict[id_row] = float(sum_row)
        race_to_id_to_sum[race] = my_dict

    # FIXED UNTIL HERE
    # race -> id_allele1_allele2 -> dict {id: allele_1: allele_2: probability:}
    race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict = {}

    with open(data_file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index % 1000 == 0:
                print(f'index: {index}, level_: {level_}')
            lst = line.split(',')
            # '47919664'
            id_row = lst[0]

            # '2.16e-13'
            probability = float(lst[2])

            level_position = LEVELS_DICT[level_]

            # A*33:01~B*14:01~C*08:02~DQB1*05:01~DRB1*01:02+A*02:01~B*39:01~C*12:03~DQB1*05:03~DRB1*14:01
            alleles_full_string = lst[1].split('+')

            left_allele = alleles_full_string[0].split('~')[level_position]
            right_allele = alleles_full_string[1].split('~')[level_position]

            # left_allele = alleles_full_string[level_position].split(',')[0]
            # right_allele = alleles_full_string[level_position].split(',')[1]

            race = _get_race(id_row)
            if race not in races_list:
                continue
            probability = float(probability) / race_to_id_to_sum[race][id_row]

            if race not in race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict:
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race] = {}

            left, right = min(left_allele, right_allele), max(left_allele, right_allele)
            # left, right = min(int(left_allele), int(right_allele)), max(int(left_allele), int(right_allele))
            id_allele1_allele2 = id_row + '_' + str(left) + '_' + str(right)

            # if it's the first time we observe this person with those alleles
            if id_allele1_allele2 not in race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race]:
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race][id_allele1_allele2] = {
                    'id': id_row,
                    'allele_1': left,
                    'allele_2': right,
                    'probability': float(probability)
                }
            else:
                # we already observed this person with those alleles
                race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race][id_allele1_allele2][
                    'probability'] += probability

    for race in races_list:
        # id_allele1_allele2 -> {id, allele_1, allele_2, probability}
        my_dict = race_to_id_allele1_allele2_to_id_allele1_allele2_probability_dict[race]
        # create a dataframe from a list of dicts
        df = pd.DataFrame.from_dict(list(my_dict.values()))
        # save the dataframe to a csv file
        df.to_csv(os.path.join(directory_path, 'levels', level_, 'races', race, 'id_allele1_allele2_probability'),
                  index=False)


def _save_results(directory_path=''):
    # directory_path = f'../data/sapir_data/{directory_name}'
    results_path = os.path.join(directory_path, "results")
    plots_path = os.path.join(directory_path, "plots")

    levels_list = LEVELS_LIST
    races_list = RACES_1_LIST

    columns = list(races_list)
    index = list(levels_list)

    df_asta_results = pd.DataFrame(columns=columns, index=index)
    df_asta = pd.DataFrame(columns=columns, index=index)
    df_dof = pd.DataFrame(columns=columns, index=index)

    for i, level in enumerate(index):
        for j, race in enumerate(columns):
            file_path = os.path.join(directory_path, 'levels', level, 'races', race, 'id_allele1_allele2_probability')
            print(f'level: {level}, race: {race}')
            p_val, statistic, dof = asta.full_algorithm(file_path=file_path,
                                                        is_first_row_contains_columns_names=True,
                                                        cutoff_value=2.0,
                                                        should_save_csv=False,
                                                        should_save_plot=os.path.join(plots_path,
                                                                                      f'level_{level}_plot.svg'),
                                                        title=f'Level {level}')
            plt.clf()
            plt.close()

            df_asta_results.iloc[i, j] = p_val
            df_asta.iloc[i, j] = statistic
            df_dof.iloc[i, j] = dof

    df_asta_results.to_csv(os.path.join(results_path, 'p_values'))
    df_asta.to_csv(os.path.join(results_path, 'statistic'))
    # df_ambiguities.to_csv(f'{real_data_path}/results/ambiguities')
    df_dof.to_csv(os.path.join(results_path, 'dof'))


def _get_race(id_row):
    return 'ALL'


def save_plots_and_data(data_file_path: str, directory_path: str = ''):
    """
    This function creates directories, sums files, probability files save the results, and finally save the plots,
    :param data_file_path: A path to your freqs file.
    The freqs file should be in the following format:
    553,A*01:01~B*49:01~C*07:01~DQB1*02:01~DRB1*07:01+A*03:01~B*07:02~C*07:02~DQB1*03:01~DRB1*11:04,8.807681966366988e-06,0
    553,A*03:01~B*49:01~C*07:01~DQB1*02:01~DRB1*07:01+A*01:01~B*07:02~C*07:02~DQB1*03:01~DRB1*11:04,2.7491553814098375e-07,1
    553,A*03:01~B*07:02~C*07:02~DQB1*02:01~DRB1*07:01+A*01:01~B*49:01~C*07:01~DQB1*03:01~DRB1*11:04,1.2764664565168216e-08,2
    207,A*01:01~B*58:01~C*07:01~DQB1*02:01~DRB1*07:01+A*02:01~B*38:01~C*12:03~DQB1*06:03~DRB1*13:01,1.682524331112655e-05,0

    Notice that the last value is the person index. The person index is 0 if it's the first time we see this person.
    Also, it is assumed that the same ids are grouped together.

    :param directory_path: A directory path (including the name of the directory) of where you want to store all the results and plots
    :return:
    """
    _create_directories(directory_path=directory_path)  # can insert path to the directory
    print("Directories created successfully")

    _create_sums_files(data_file_path, directory_path=directory_path)
    print("Sums files created successfully")

    # create for each level the p_k(i,j) files
    for level in LEVELS_LIST:
        _create_id_allele1_allele2_probability_files_for_each_level(level,
                                                                    data_file_path,
                                                                    directory_path=directory_path)
    print("Probability files created successfully")

    # save results for each level
    _save_results(directory_path=directory_path)
    print("Results saved successfully")
