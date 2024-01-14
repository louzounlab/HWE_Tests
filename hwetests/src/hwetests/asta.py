import numpy as np
import math
import pandas as pd
import scipy.stats as stats
import csv
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import RATIO


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(alleles_amount, population_amount_, alleles_probabilities,
                                observed_probabilities, correction, cutoff):
    value = 0.0
    amount_of_small_expected_ = 0
    for row in range(alleles_amount):
        for col in range(row, alleles_amount):
            mult = 1.0
            if row != col:
                mult = 2.0
            expected_val = mult * population_amount_ * alleles_probabilities[row] * alleles_probabilities[col]
            observed_val = population_amount_ * observed_probabilities[row, col]
            correction_val = correction[row, col]
            variance_val = expected_val * correction_val
            if expected_val < cutoff:
                amount_of_small_expected_ += 1
                continue
            value += ((expected_val - observed_val) ** 2) / variance_val
    return value, amount_of_small_expected_


def run_experiment(alleles_count, population_amount, alleles_probabilities,
                   observed_probabilities, correction, var_obs, index_to_allele_, should_save_csv_, cutoff_value_):
    chi_squared_stat, amount_of_small_expected = calculate_chi_squared_value(alleles_amount=alleles_count,
                                                                             population_amount_=population_amount,
                                                                             alleles_probabilities=alleles_probabilities,
                                                                             observed_probabilities=observed_probabilities,
                                                                             correction=correction,
                                                                             cutoff=cutoff_value_)
    couples_amount = (alleles_count * (alleles_count + 1)) / 2 - 1
    dof = couples_amount - amount_of_small_expected

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)

    if should_save_csv_:
        if isinstance(should_save_csv_, str):
            file_name = should_save_csv_ + '.csv'
        else:
            file_name = 'alleles_data.csv'
        columns = ['first_allele', 'second_allele', 'O(i,j)', 'E(i,j)', 'Rho(i,j)', 'E(i,j)*Rho(i,j)',
                   'var_observed(i,j)']
        with open(file_name, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(columns)

            for i in range(alleles_count):
                for j in range(i, alleles_count):
                    mult = 1.0
                    if i != j:
                        mult = 2.0
                    expected_val = mult * population_amount * alleles_probabilities[i] * alleles_probabilities[j]
                    observed_val = population_amount * observed_probabilities[i, j]
                    correction_val = correction[i, j]
                    variance_val = expected_val * correction_val
                    var_observed_value = var_obs[i, j]

                    first_allele = index_to_allele_[i]
                    second_allele = index_to_allele_[j]

                    writer.writerow([first_allele, second_allele, observed_val, expected_val, correction_val,
                                     variance_val,
                                     var_observed_value])

    return p_value, chi_squared_stat, dof


def full_algorithm(file_path,
                   is_first_row_contains_columns_names=False,
                   cutoff_value=0.0,
                   should_save_csv=False,
                   should_save_plot=False,
                   title=''):
    """
    ASTA Algorithm.

    Performs a modified Chi-Squared statistical test on ambiguous observations.
    :param file_path: A path to a csv file with columns: 1) index or id of a donor (integer or string).
    Assuming columns are separated with , or + and no whitespaces in the csv file. 2) first allele (integer or string).
    3) second allele (integer or string). 4) probability (float).
    :param is_first_row_contains_columns_names: If True then it is assumed that the first row contains the
    columns names: i.e. 'column1,column2,...'.
    :param cutoff_value: (optional, default value is 0.0) A float value that decides to not account (O-E)^2 / E in the summation
    of the Chi-Squared statistic if E < cutoff.
    :param should_save_csv: (optional, default value is False) Either boolean or string, if it's True then a csv with the columns:
     [first allele, second allele, observed, expected, variance] is saved (named 'alleles_data.csv')
    and if it's a string then a csv with the given string name is saved.
    :param should_save_plot: (optional, default value is False) Either boolean or string, if it's True then
    a png containing 2 bar plots is saved (named 'alleles_barplot.png') for each allele showing its chi squared statistic over degrees of freedom
    (summing over the observations only associated with this allele) and -log_10(p_value).
    If it's a string then a csv with the given string name is saved.
    :param title: (optional, default value is '') A string that will be the title of the plot.
    :return: p-value (float), Chi-Squared statistic (float), degrees of freedom (integer). also saves a csv.
    with columns: first allele, second allele, observed, variance:
    """
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}

    # id -> {'i_j' -> [i, j, O_k(i,j)], ...}
    # actual id string, i and j are indices
    id_to_i_j_to_i_j_observation = {}

    # first read all the rows and get indices of ids and alleles and amounts
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if (index == 0) and is_first_row_contains_columns_names:
                continue
            lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()

            id_row = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            observation = float(lst[3])

            if id_row not in id_to_index:
                id_to_index[id_row] = len(id_to_index)

            if allele_1 not in allele_to_index:
                # updating the inverse dictionary
                index_to_allele[len(allele_to_index)] = allele_1
                # updating the dictionary
                allele_to_index[allele_1] = len(allele_to_index)
            if allele_2 not in allele_to_index:
                # updating the inverse dictionary
                index_to_allele[len(allele_to_index)] = allele_2
                # updating the dictionary
                allele_to_index[allele_2] = len(allele_to_index)

            # get indices of alleles and make i<=j
            i = allele_to_index[allele_1]
            j = allele_to_index[allele_2]
            i, j = min(i, j), max(i, j)
            i_j = f'{i}_{j}'
            # update dictionary
            if id_row not in id_to_i_j_to_i_j_observation:
                id_to_i_j_to_i_j_observation[id_row] = {i_j: [i, j, observation]}
            else:
                # id exists in dictionary, we need to append the observation for i, j
                # if i_j exists, we need to add the probability and if it doesn't,
                # we need to append a new observation

                # if observation i, j doesn't exist
                if i_j not in id_to_i_j_to_i_j_observation[id_row]:
                    id_to_i_j_to_i_j_observation[id_row][i_j] = [i, j, observation]
                else:
                    # observation i, j exists. add the new observation
                    id_to_i_j_to_i_j_observation[id_row][i_j][2] += observation

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    #
    # # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # go over the ids
    for current_id in id_to_i_j_to_i_j_observation:
        # sum observations
        sum_observations = 0.0
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            sum_observations += id_to_i_j_to_i_j_observation[current_id][i_j][2]
        # get the probabilities
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            # i and j are already i <= j
            i = id_to_i_j_to_i_j_observation[current_id][i_j][0]
            j = id_to_i_j_to_i_j_observation[current_id][i_j][1]
            probability = id_to_i_j_to_i_j_observation[current_id][i_j][2] / sum_observations

            alleles_probabilities[i] += 0.5 * probability
            alleles_probabilities[j] += 0.5 * probability
            observed_probabilities[i, j] += probability
            correction[i, j] += (probability ** 2)

    # p(i) = sum_k_j p_k(i,j) / N
    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(i, alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    var_obs = np.zeros(shape=(alleles_count, alleles_count))
    if should_save_csv:
        for i in range(alleles_count):
            for j in range(i, alleles_count):
                arr = []
                for current_id in id_to_i_j_to_i_j_observation:
                    i_j = f'{i}_{j}'
                    if i_j in id_to_i_j_to_i_j_observation[current_id]:
                        arr.append(id_to_i_j_to_i_j_observation[current_id][i_j][2])
                    else:
                        arr.append(0.0)
                # calculate observed variance
                mean_value = np.mean(arr)
                var_value = 0.0
                for a in arr:
                    var_value += (a - mean_value) ** 2
                var_obs[i, j] = var_value

    p_value, chi_squared, dof = run_experiment(alleles_count=alleles_count,
                                               population_amount=population_amount,
                                               alleles_probabilities=alleles_probabilities,
                                               observed_probabilities=observed_probabilities,
                                               correction=correction,
                                               var_obs=var_obs,
                                               index_to_allele_=index_to_allele,
                                               should_save_csv_=should_save_csv,
                                               cutoff_value_=cutoff_value)

    if should_save_plot:
        # save a bar plot showing for each allele its deviation from HWE
        # couples_amount = int((alleles_count * (alleles_count + 1)) / 2 - 1)
        df = pd.DataFrame(index=range(alleles_count), columns=['Alleles', 'Normalized statistic', '-log_10(p_value)'])
        logs_list = ['s' for _ in range(alleles_count)]
        p_values = [0.0 for _ in range(alleles_count)]
        for i in range(alleles_count):
            # for allele i: calculate Statistic and p_value
            statistic_i = 0.0
            amount_of_small_expected_i = 0
            for j in range(alleles_count):
                t, m = min(i, j), max(i, j)
                mult = 1
                if t != m:
                    mult = 2
                expected_val_i = mult * population_amount * alleles_probabilities[t] * alleles_probabilities[m]
                observed_val_i = population_amount * observed_probabilities[t, m]
                correction_val_i = correction[t, m]
                variance_val_i = expected_val_i * correction_val_i
                if variance_val_i < cutoff_value:
                    amount_of_small_expected_i += 1
                    continue
                statistic_i += ((expected_val_i - observed_val_i) ** 2) / variance_val_i
            # calculate degrees of freedom
            dof_i = (alleles_count - 1) - amount_of_small_expected_i

            # if degrees of freedom for allele i is too small and we have many alleles, don't plot this allele
            if (dof_i < 10) and (alleles_count > 11):
                df.loc[i] = ['', 0, 's']
                continue

            allele_i = index_to_allele[i]

            p_value_i = 1 - stats.chi2.cdf(x=statistic_i,
                                           df=dof_i)
            p_values.append(p_value_i)

            if p_value_i == 0.0:
                logs_list[i] = 50
            else:
                logs_list[i] = -math.log(p_value_i, 10)
            df.loc[i] = [allele_i, statistic_i / dof_i, logs_list[i]]
        # we might have a -log(p_value) bigger than 50, therefore we need to find the max and update the df
        logs_list = [x for x in logs_list if (type(x) == int or type(x) == float)]
        if not logs_list:
            print('All the alleles have a small dof while having many alleles, therefore the plot will not be shown')
            return p_value, chi_squared, dof
        max_log = max(logs_list)
        if max_log > 50:
            for i in range(alleles_count):
                if p_values[i] == 0.0:
                    df.loc[i, '-log_10(p_value)'] = max_log

        # sort dataframe according to p_values and take the smallest 20 statistics.
        df = df.loc[df['-log_10(p_value)'] != 's']
        df = df.sort_values('Normalized statistic', ascending=False).head(min(alleles_count, 20))
        # we need to update the log p-values (some may be infinite, so we set them to the max value from the 20)

        # plot the dataframe into 2 bar plots
        # fig, axes = plt.subplots(1, 2)
        # plt.subplot(1, 2, 1)
        # sns.set_color_codes('pastel')
        # sns.barplot(x='-log_10(p_value)', y='Alleles', data=df,
        #             label='Total', color='royalblue', edgecolor='w')
        # sns.set_color_codes('muted')
        # # invert the order of x-axis values
        # ax = plt.gca()
        # ax.set_xlim(ax.get_xlim()[::-1])
        # # ax.yaxis.set_label_position("right")
        # ax.yaxis.tick_right()
        # plt.ylabel('')

        # move ticks to the right

        # plt.subplot(1, 2, 2)

        # making a scatter plot with color bar
        sns.set_style('white')
        plt.rcParams["font.family"] = "Times New Roman"

        # plt.figure(figsize=((6, 6)))
        plot = plt.scatter(df['Normalized statistic'], df['Alleles'], c=df['-log_10(p_value)'], cmap='Reds',
                           vmin=0, vmax=max_log)
        # plotting the color bar
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label('-log_10(p_value)', fontsize=16)
        # removing the points from the color bar
        plot.remove()
        if df['-log_10(p_value)'].nunique() == 1:
            ax = sns.barplot(x='Normalized statistic', y='Alleles', data=df,
                             hue=df['-log_10(p_value)'], palette=sns.color_palette(['#6d010e']), dodge=False)
        else:
            ax = sns.barplot(x='Normalized statistic', y='Alleles', data=df,
                             hue=df['-log_10(p_value)'], palette='Reds', dodge=False)
        # plt.legend(ncol=2, loc='lower right')
        sns.despine(left=True, bottom=True)
        ax.legend_.remove()
        # fig.tight_layout()

        # take care of the font sizes
        ax.set_xlabel('Normalized statistic', fontsize=18)
        ax.set_ylabel('Alleles', fontsize=18)
        ax.tick_params(labelsize=16)
        if isinstance(should_save_plot, str):
            file_name = should_save_plot
        else:
            file_name = 'alleles_barplot.pdf'
        if title:
            plt.title(title, fontsize=20)
        plt.savefig(file_name, format='pdf', bbox_inches="tight")
        # clear plt
        # plt.clf()
        # plt.close()

    return p_value, chi_squared, dof
