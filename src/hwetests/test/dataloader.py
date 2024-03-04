import os


def get_path(filename):
    '''
    This function returns theabsolute  path of the file
    :param filename: the name of the file
    :return: the absolute path of the file
    '''
    # the absolute path of this file
    current_path = os.path.dirname(__file__)
    final_path = os.path.join(current_path, "data", filename)
    return final_path