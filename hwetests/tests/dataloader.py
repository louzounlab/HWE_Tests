import os


def get_path(is_ambiguous: bool):
    """
    This function returns the absolute path of the file :param is_ambiguous: a boolean value to indicate if we should
    return the path of the ambiguous or unambiguous data file.
    :return: the absolute path of the file
    """
    if is_ambiguous:
        filename = "ambiguous_data.csv"
    else:
        filename = "unambiguous_data.csv"
    # the absolute path of this file
    current_path = os.path.dirname(__file__)
    final_path = os.path.join(current_path, "data", filename)
    return final_path
