import numpy as np
import math
import pandas as pd
import scipy.stats as stats
import csv
from matplotlib import pyplot as plt
import seaborn as sns


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
            if variance_val < cutoff:
                amount_of_small_expected_ += 1
                continue
            value += ((expected_val - observed_val) ** 2) / variance_val
    return value, amount_of_small_expected_


def run_experiment(alleles_count, population_amount, alleles_probabilities,
                   observed_probabilities, correction, index_to_allele_, should_save_csv_, cutoff_value_):
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
        columns = ['first_allele', 'second_allele', 'observed', 'expected', 'variance']
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

                    first_allele = index_to_allele_[i]
                    second_allele = index_to_allele_[j]

                    writer.writerow([first_allele, second_allele, observed_val, expected_val, variance_val])

    return p_value, chi_squared_stat, dof


def full_algorithm(file_path, cutoff_value=0.0, should_save_csv=False, should_save_plot=False):
    """
    ASTA Algorithm.

    Performs a modified Chi-Squared statistical test on ambiguous observations.
    :param file_path: A path to a csv file with columns: 1) index or id of a donor (integer or string).
    2) first allele (integer or string). 3) second allele (integer or string). 4) probability (float).
    :param cutoff_value: (optional, default value is 0.0) A float value that decides to not account (O-E)^2 / E in the summation
    of the Chi-Squared statistic if E < cutoff.
    :param should_save_csv: (optional, default value is False) Either boolean or string, if it's True then a csv with the columns:
     [first allele, second allele, observed, expected, variance] is saved (named 'alleles_data.csv')
    and if it's a string then a csv with the given string name is saved.
    :param should_save_plot: (optional, default value is False) Either boolean or string, if it's True then
    a png containing 2 bar plots is saved (named 'alleles_barplot.png') for each allele showing its chi squared statistic over degrees of freedom
    (summing over the observations only associated with this allele) and -log_10(p_value).
    If it's a string then a csv with the given string name is saved.
    :return: p-value (float), Chi-Squared statistic (float), degrees of freedom (integer). also saves a csv.
    with columns: first allele, second allele, observed, variance:
    """
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}
    # index1_index2 -> [obs_prob, corr]

    # first read all the rows and get indices of ids and alleles and amounts
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            # probability = float(lst[3])

            if id not in id_to_index:
                id_to_index[id] = len(id_to_index)

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

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))
    #
    # # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # calculate {p_k(i,j)}
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            # id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            probability = float(lst[3])

            # id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)

            alleles_probabilities[allele_1_index] += 0.5 * probability
            alleles_probabilities[allele_2_index] += 0.5 * probability

            observed_probabilities[allele_1_index, allele_2_index] += probability

            correction[allele_1_index, allele_2_index] += (probability ** 2)

    # p(i) = sum_k_j p_k(i,j) / N
    alleles_probabilities /= population_amount

    for i in range(alleles_count):
        for j in range(alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    # for i in range(alleles_count):
    #     for j in range(alleles_count):
    #         observed_probabilities[i, j] /= population_amount
    #         if observed_probabilities[i, j] == 0:
    #             correction[i, j] = 1.0
    #         else:
    #             correction[i, j] /= (population_amount * observed_probabilities[i, j])

    p_value, chi_squared, dof = run_experiment(alleles_count=alleles_count,
                                               population_amount=population_amount,
                                               alleles_probabilities=alleles_probabilities,
                                               observed_probabilities=observed_probabilities,
                                               correction=correction,
                                               index_to_allele_=index_to_allele,
                                               should_save_csv_=should_save_csv,
                                               cutoff_value_=cutoff_value)

    if should_save_plot:
        # save a bar plot showing for each allele its deviation from HWE
        couples_amount = int((alleles_count * (alleles_count + 1)) / 2 - 1)
        df = pd.DataFrame(index=range(alleles_count), columns=['Alleles', 'Normalized statistic', '-log_10(p_value)'])
        logs_list = [0 for _ in range(alleles_count)]
        for i in range(alleles_count):
            # for allele i: calculate Statistic and p_value
            statistic = 0.0
            amount_of_small_expected = 0
            for j in range(alleles_count):
                t, m = min(i, j), max(i, j)
                mult = 1
                if t != m:
                    mult = 2
                expected_val = mult * population_amount * alleles_probabilities[t] * alleles_probabilities[m]
                observed_val = population_amount * observed_probabilities[t, m]
                correction_val = correction[t, m]
                variance_val = expected_val * correction_val
                if variance_val < cutoff_value:
                    amount_of_small_expected += 1
                    continue
                statistic += ((expected_val - observed_val) ** 2) / variance_val
            # calculate degrees of freedom
            dof_i = (alleles_count - 1) - amount_of_small_expected

            allele = index_to_allele[i]

            if dof_i <= 0:
                df.loc[i] = [allele, 2.0, '?']
                continue
            # calculate p_value
            p_value_i = 1 - stats.chi2.cdf(x=statistic,
                                           df=dof_i)
            if p_value_i == 0.0:
                logs_list[i] = 'infty'
            else:
                logs_list[i] = -math.log(p_value_i, 10)
            df.loc[i] = [allele, statistic / dof, logs_list[i]]

        # sort dataframe according to p_values and take the smallest 20 statistics.
        df = df.sort_values('Normalized statistic').head(min(alleles_count, 20))
        # we need to update the log p-values (some may be infinite, so we set them to the max value from the 20)
        logs_list_ints = [logs_list[i] for i in df.index if isinstance(logs_list[i], (int, float))]
        max_log = max(logs_list_ints)
        for i in df.index:
            if logs_list[i] == 'infty':
                logs_list[i] = max_log
            df.loc[i, '-log_10(p_value)'] = logs_list[i]
        df = df.loc[df['-log_10(p_value)'] != '?']
        # plot the dataframe into 2 bar plots
        fig, axes = plt.subplots(1, 2)
        plt.subplot(1, 2, 1)
        sns.set_color_codes('pastel')
        sns.barplot(x='-log_10(p_value)', y='Alleles', data=df,
                    label='Total', color='royalblue', edgecolor='w')
        sns.set_color_codes('muted')
        # invert the order of x-axis values
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
        # ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('')

        # move ticks to the right

        plt.subplot(1, 2, 2)
        sns.barplot(x='Normalized statistic', y='Alleles', data=df,
                    color='slategrey', edgecolor='w')
        # plt.legend(ncol=2, loc='lower right')
        sns.despine(left=True, bottom=True)
        fig.tight_layout()

        if isinstance(should_save_plot, str):
            file_name = should_save_plot + '.png'
        else:
            file_name = 'alleles_barplot.png'
        plt.savefig(file_name, pad_inches=0.2, bbox_inches="tight")

    return p_value, chi_squared, dof
