import numpy as np
import pandas as pd
import scipy.stats as stats
import random


# need to calculate only on the upper triangle because the matrices are symmetric
def calculate_chi_squared_value(population_amount_, counts_expected_, observed_probabilities_, correction_):
    value = 0
    for row in range(counts_expected_.shape[0]):
        for col in range(row, counts_expected_.shape[1]):
            expected_ = counts_expected_[row, col]
            observed_ = population_amount_ * observed_probabilities_[row, col]
            variance_ = expected_
            # value += correction_[row, col] * (((expected_ - observed_) ** 2) / variance_)
            value += (1 / correction_[row, col]) * (((expected_ - observed_) ** 2) / variance_)
    return value


def run_experiment(alleles_count, population_amount, alleles_probabilities,
                   observed_probabilities, correction):
    counts_expected = np.zeros((alleles_count, alleles_count))
    for j in range(alleles_count):
        for k in range(j, alleles_count):
            # EVEN
            mult = 1
            if k != j:
                mult = 2
            expected_value = mult * population_amount * alleles_probabilities[k] * alleles_probabilities[j]
            counts_expected[k, j] = expected_value
            counts_expected[j, k] = expected_value

    chi_squared_stat = calculate_chi_squared_value(population_amount_=population_amount,
                                                   counts_expected_=counts_expected,
                                                   observed_probabilities_=observed_probabilities,
                                                   correction_=correction)
    dof = (alleles_count * (alleles_count + 1)) / 2 - 1

    # print(f' alpha for choice: {alpha_val}')
    # print(f' chi square value: {chi_squared_stat}')

    # crit = stats.chi2.ppf(q=0.95, df=dof)
    # print(f'Critical value: {crit}')

    p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,
                                 df=dof)
    return p_value, chi_squared_stat, dof


def full_algorithm(file_path):
    """
        Modified Chi-Squared Algorithm.

        Performs a modified Chi-Squared statistical test on ambiguous observations.
        :param file_path: A path to a csv file with columns: 1) index or id of a donor.
        2) first allele (integer). 3) second allele (integer). 4) probability (float).
        :return: A p-value
        """
    id_to_index = {}
    allele_to_index = {}

    # first read all the rows and get indices of ids and alleles and amounts
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            id = lst[0]
            allele_1 = int(lst[1])
            allele_2 = int(lst[2])
            allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            probability = float(lst[3])

            if id not in id_to_index:
                id_to_index[id] = len(id_to_index)

            if allele_1 not in allele_to_index:
                allele_to_index[allele_1] = len(allele_to_index)
            if allele_2 not in allele_to_index:
                allele_to_index[allele_2] = len(allele_to_index)

    alleles_count = len(allele_to_index)
    population_amount = len(id_to_index)

    # {p(i)}
    alleles_probabilities = np.zeros(alleles_count)

    # {p(i, j)}
    observed_probabilities = np.zeros(shape=(alleles_count, alleles_count))

    # correction matrix
    correction = np.zeros(shape=(alleles_count, alleles_count))

    # calculate {p_k(i,j)}
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            id = lst[0]
            allele_1 = int(lst[1])
            allele_2 = int(lst[2])
            allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            probability = float(lst[3])

            id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            alleles_probabilities[allele_1_index] += 0.5 * probability
            alleles_probabilities[allele_2_index] += 0.5 * probability

            observed_probabilities[allele_1_index, allele_2_index] += probability

            correction[allele_1_index, allele_2_index] += (probability ** 2)

    for i in range(alleles_count):
        for j in range(alleles_count):
            observed_probabilities[i, j] /= population_amount
            if observed_probabilities[i, j] == 0:
                correction[i, j] = 1.0
            else:
                correction[i, j] /= (population_amount * observed_probabilities[i, j])

    p_value, chi_squared, dof = run_experiment(alleles_count=alleles_count,
                                               population_amount=population_amount,
                                               alleles_probabilities=alleles_probabilities,
                                               observed_probabilities=observed_probabilities,
                                               correction=correction)

    return p_value
