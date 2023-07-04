"""
General functions used to read YAML files
"""

# pylint: disable=import-error

import yaml

class loadYAML():
    """
    Object which represents the contents of a YAML file in
    numpy objects which can be more easily manipulated in the
    rest of the code

    Parameters
    ----------
    filenames : str/list of str
                The name of the YAML file to be represented
    """

    def __init__(self, filename):

        with open(filename, "r", encoding="utf-8") as stream:
            input_data = yaml.safe_load(stream)

        self.keys = []
        for key_name in input_data.keys():
            self.keys.append(key_name)
            setattr(self, key_name, input_data[key_name])
    #end def __init__
#end class loadYAML
