import numpy as np
import pandas as pd
import random
import umat


def full_algorithm(file_path):
    """
    UMAT with sampling algorithm.

    Performs UMAT test on ambiguous observations.
    :param file_path: A path to a csv file with columns: 1) index or id of a donor (integer or string).
    2) first allele (integer or string). 3) second allele (integer or string). 4) probability (float).+
    :return: p-value (float), degrees of freedom (integer), Chi-Squared statistic (float), also saves a csv.
    with columns: first allele, second allele, observed, variance:
    """
    id_to_index = {}
    allele_to_index = {}
    index_to_allele = {}

    # first read all the rows and get indices of ids and alleles and amounts
    with open(file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            # header
            if index == 0:
                columns = line.strip('\n').split(',')
                id_col = columns[0]
                allele_1_col = columns[1]
                allele_2_col = columns[2]
                probability_col = columns[3]
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

    # {O(i, j)}
    observations = np.zeros(shape=(alleles_count, alleles_count))

    # sort by ids
    df = pd.read_csv(file_path,
                     dtype={id_col: str,
                            allele_1_col: str,
                            allele_2_col: str,
                            probability_col: float})
    df = df.sort_values(by=id_col)
    temp_file_path = 'temp_file.csv'
    df.to_csv(temp_file_path)

    # now the ids are sorted, for each donor k we sample alleles (i, j) from his set of probabilities
    # the last id that we scanned
    last_id = ''
    # the probabilities of the last id. list of floats
    last_probabilities = []
    # the possible alleles of the last id, according to his probabilities. list where each element is a list [i, j]
    last_possible_alleles = []
    with open(temp_file_path, encoding="utf8") as infile:
        for index, line in enumerate(infile):
            if index == 0:
                continue
            lst = line.strip('\n').split(',')

            current_id = lst[0]

            # if we finished scanning the last id
            if last_id and (current_id != last_id):
                # sample alleles for the last id
                index = random.choices(population=range(len(last_probabilities)), weights=last_probabilities,
                                       k=1)[0]
                # get the sampled alleles
                i, j = last_possible_alleles[index]
                # update observations
                observations[i, j] += 1
                # reset last person
                last_id = current_id
                last_probabilities = []
                last_possible_alleles = []

            allele_1 = lst[1]
            allele_2 = lst[2]
            # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            probability = float(lst[3])

            id_index = id_to_index[id]

            allele_1_index = allele_to_index[allele_1]
            allele_2_index = allele_to_index[allele_2]

            allele_1_index, allele_2_index = min(allele_1_index, allele_2_index), max(allele_1_index, allele_2_index)

            last_probabilities.append(probability)
            last_possible_alleles.append([allele_1_index, allele_2_index])

    # we still have the last id
    index = random.choices(population=range(len(last_probabilities)), weights=last_probabilities,
                           k=1)[0]
    # get the sampled alleles
    i, j = last_possible_alleles[index]
    # update observations
    observations[i, j] += 1
    # reset last person
    last_id = current_id
    last_probabilities = []
    last_possible_alleles = []

    return umat.full_algorithm(observations=observations)
