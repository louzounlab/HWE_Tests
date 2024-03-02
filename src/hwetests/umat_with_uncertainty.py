import numpy as np
import random
from . import umat


def full_algorithm(file_path, is_first_row_contains_columns_names=False,
                   start_from=None, iterations=None):
    """
    UMAT with sampling algorithm.

    Performs UMAT test on ambiguous observations.
    :param file_path: A path to a csv file with columns: 1) index or id of a donor (integer or string).
    2) first allele (integer or string). 3) second allele (integer or string). 4) probability (float).
    Assuming columns are separated with , or + and no whitespaces in the csv file
    :param is_first_row_contains_columns_names: If True then it is assumed that the first row contains the
    columns names: i.e. 'column1,column2,...'
    :param start_from:
    :param iterations:
    :return: A p-value under the null Hypothesis that observations are distributed around HWE.
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
                # columns = line.strip('\n').split(',')
                # id_col = columns[0]
                # allele_1_col = columns[1]
                # allele_2_col = columns[2]
                # probability_col = columns[3]
                continue
            lst = line.strip('\n').replace('+', ' ').replace(',', ' ').split()

            current_id = lst[0]
            allele_1 = lst[1]
            allele_2 = lst[2]
            # allele_1, allele_2 = min(allele_1, allele_2), max(allele_1, allele_2)
            observation = float(lst[3])

            if current_id not in id_to_index:
                id_to_index[current_id] = len(id_to_index)

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
            if current_id not in id_to_i_j_to_i_j_observation:
                id_to_i_j_to_i_j_observation[current_id] = {i_j: [i, j, observation]}
            else:
                # id exists in dictionary, we need to append the observation for i, j
                # if i_j exists, we need to add the probability and if it doesn't,
                # we need to append a new observation

                # if observation i, j doesn't exist
                if i_j not in id_to_i_j_to_i_j_observation[current_id]:
                    id_to_i_j_to_i_j_observation[current_id][i_j] = [i, j, observation]
                else:
                    # observation i, j exists. add the new observation
                    id_to_i_j_to_i_j_observation[current_id][i_j][2] += observation

    alleles_count = len(allele_to_index)
    # population_amount = len(id_to_index)

    # {O(i, j)}
    observations = np.zeros(shape=(alleles_count, alleles_count))

    # go over the ids
    for current_id in id_to_i_j_to_i_j_observation:
        # list of observations of current id
        observed_list = []
        # list if indices i_j of current id
        indices_list = []
        # sum of observations of current id
        sum_observed = 0.0
        for i_j in id_to_i_j_to_i_j_observation[current_id]:
            indices_list.append(i_j)
            observed_list.append(id_to_i_j_to_i_j_observation[current_id][i_j][2])
            sum_observed += id_to_i_j_to_i_j_observation[current_id][i_j][2]
        # normalize observations of current id into probabilities
        for t in range(len(observed_list)):
            observed_list[t] /= sum_observed
        # sample from probabilities of current id and get the corresponding alleles i, j
        index = random.choices(population=range(len(observed_list)), weights=observed_list,
                               k=1)[0]
        i_j = indices_list[index]
        i = id_to_i_j_to_i_j_observation[current_id][i_j][0]
        j = id_to_i_j_to_i_j_observation[current_id][i_j][1]

        # update observations
        observations[i, j] += 1

    return umat.full_algorithm(observations=observations,
                               start_from=start_from,
                               iterations=iterations)
