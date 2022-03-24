# SPDX-License-Identifier: MIT

""" Handle interactions with HDF5
"""

import h5py as h5
import numpy as np

class HDF5Group():
    """Object which represents the contents of a HDF5 group"""
    def __init__(self):
        self.keys = []

    def add_attr(self, name, val):
        """ Add a new attribute to the group
        """
        self.keys.append(name)
        setattr(self, name, val)

#--------------------------------------------------
class loadHDF5():
    """Object which represents the contents of a HDF5 file in
    numpy objects which can be more easily manipulated in the
    rest of the code

    Parameters
    ----------
    filenames : str/list of str
                The name of the HDF5 files to be represented
    name_map  : dict
                A dictionary indicating any desired name remapping
    """

    def __init__(self, filenames, name_map = None):
        if isinstance(filenames, list):
            nb_diag = len(filenames)
        else:
            nb_diag = 1
            filenames = [filenames]

        if name_map is None:
            name_map = dict()

        self.keys = []

        for i,f in enumerate(filenames):

            print('---> Reading of ' + str(f))
            with h5.File(f, 'r') as fh5:

                if nb_diag == 1:
                    self.collect_group(fh5, name_map, self)
                else:
                    self.collect_time_group(fh5, name_map, i, nb_diag, self)
            #end with fh5
        #end for
    #end def __init__

    @staticmethod
    def get_str_var_new(str_var, name_map):
        """ replace '%' by '_' in variable names

        Parameters
        ----------
        str_var  : str
                   The variable name
        name_map : dict
                   A dictionary mapping variable names to their replacement
        """
        str_var_new = str_var.replace('%','_')
        if str_var_new in name_map:
            str_var_new = name_map[str_var_new]
        return str_var_new

    def collect_group(self, fh5, name_map, group):
        """ Collect all elements from a group from a HDF5 file
        """
        for str_var in list(fh5.keys()):
            str_var_new = self.get_str_var_new(str_var, name_map)
            var_tmp     = fh5[str_var]

            if isinstance(var_tmp, h5.Group):
                var = HDF5Group()
                # Collect subgroup
                self.collect_group(var_tmp, name_map, var)
            else:
                # Collect data
                if var_tmp.shape == ():
                    var = var_tmp[()]
                    if isinstance(var, np.bytes_):
                        var = var.decode('UTF-8')
                elif var_tmp.size == 1:
                    var = var_tmp[0]
                    if isinstance(var, np.generic):
                        var = var.item()
                else:
                    var = np.array(var_tmp, order='C')

            group.add_attr(str_var_new, var)

    def collect_time_group(self, fh5, name_map, i_diag, nb_diag, group):
        """ Collect all elements from a group from a HDF5 file.
        This function should be called multiple times to create an array of
        results from files representing different points in time

        Parameters
        ----------
        fh5      : h5.File
                   The HDF5 file to be read
        name_map : dict
                   A dictionary indicating if and and how any attributes
                   should be renamed
        i_diag   : int
                   The index of the current diagnostic file
        nb_diag  : int
                   The total number of current diagnostic files
        group    : HDF5Group
                   A group where the results will be stored
        """

        for str_var in list(fh5.keys()):
            str_var_new = self.get_str_var_new(str_var, name_map)
            var_tmp     = fh5[str_var]

            nd = len(var_tmp.shape)

            if i_diag == 0:
                # Create enough space for data at each time step
                if isinstance(var_tmp, h5.Group):
                    var = HDF5Group()
                else:
                    # Fortran ordering is used as we index by time a lot

                    var = np.empty(shape = (nb_diag,*var_tmp.shape),
                                   dtype = var_tmp.dtype,
                                   order='C')
                group.add_attr(str_var_new, var)
            else:
                # Collect previously created space
                var = getattr(group, str_var_new)

            if isinstance(var_tmp, h5.Group):
                # Collect subgroup
                self.collect_time_group(var_tmp, name_map, i_diag, nb_diag, var)
            else:
                # Create indexes, last index is time
                idxs = [i_diag] + [slice(None)]*nd
                # Save data
                var[tuple(idxs)] = var_tmp[()]
        #end for

    def append(self, other_file):
        """ Save the values from another loadHDF5 instance into this instance
        """
        for k in other_file.keys:
            self.add_attr(k, getattr(other_file,k))

    def add_attr(self, name, val):
        """ Add a new attribute to the file
        """
        self.keys.append(name)
        setattr(self, name, val)

#end class loadHDF5

#--------------------------------------------------
def write_group(myfile, group):
    """ Save the group to a HDF5 file
    """
    for name in group.keys:
        data = getattr(group, name)

        if isinstance(data, HDF5Group):
            subgroup = myfile.create_group(name)
            write_group(subgroup, data)
        else:

            if hasattr(data, 'shape'):
                dset = myfile.create_dataset(name, data.shape, dtype =data.dtype )
                dset[()] = data
            else:
                dtype = h5.string_dtype(encoding='ascii') if isinstance(data, str) else type(data)
                dset = myfile.create_dataset(name, [1], dtype = dtype, data = data)


#--------------------------------------------------
def write_results(filename, data):
    """ Save data to a HDF5 file
    """

    myfile = h5.File(filename,'w')
    write_group(myfile, data)
    myfile.close()


#-----------------------------------------------
def save_dict_contents_to_group( h5file, path, dic):
    """ Save a dictionary into an HDF5 group
    """

    # argument type checking
    if not isinstance(dic, dict):
        raise ValueError("must provide a dictionary")

    if not isinstance(path, str):
        raise ValueError("path must be a string")

    if not isinstance(h5file, h5.File):
        raise ValueError("must be an open hdf5 file")

    # save items to the hdf5 file
    for key, item in dic.items():
        #print(key,item)
        key = str(key)
        if isinstance(item, list):
            item = np.array(item)
            #print(item)
        if not isinstance(key, str):
            raise ValueError("dict keys must be strings to save to hdf5")
        # save strings, numpy.int64, and numpy.float64 types
        if isinstance(item, (np.int64, np.float64, str, np.float, float, np.float32,int)):
            #print( 'here' )
            h5file[path + key] = item
            if not h5file[path + key][()] == item:
                raise ValueError('The data representation in the HDF5 file does not match the original dict.')
        # save numpy arrays
        elif isinstance(item, np.ndarray):
            try:
                h5file[path + key] = item
            except KeyError:
                item = np.array(item).astype('|S9')
                h5file[path + key] = item
            if not np.array_equal(h5file[path + key][()], item):
                raise ValueError('The data representation in the HDF5 file does not match the original dict.')
        # save dictionaries
        elif isinstance(item, dict):
            save_dict_contents_to_group(h5file, path + key + '/', item)
        # other types cannot be saved and will result in an error
        else:
            raise ValueError('Cannot save %s type.' % type(item))

#end def save_dict_contents_to_group
#-----------------------------------------------


#-----------------------------------------------
def load_dict_contents_from_group(h5file, path):
    """ Load a dictionary from an HDF5 group
    """

    ans = {}
    for key,item in h5file[path].items():
        if isinstance(item, h5.Dataset):
            ans[key] = item[()]
        elif isinstance(item, h5.Group):
            ans[key] = load_dict_contents_from_group(h5file, path + key + '/')
    return ans

#end load_dict_contents_from_group
#-----------------------------------------------


#-----------------------------------------------
def save_dict_to_hdf5(dic, group, filename):
    """ Save a dictionary to a HDF5 file
    """

    with h5.File(filename, 'a') as h5file:
        if group not in h5file.keys():
            save_dict_contents_to_group(h5file, group, dic)
        else:
            print(group + ' already saved in ' + filename)

#end def save_dict_to_hdf5
#-----------------------------------------------


#-----------------------------------------------
def load_dict_from_hdf5(filename, group):
    """ Load a dictionary into a HDF5 file
    """

    with h5.File(filename, 'r') as h5file:
        if group in h5file.keys():
            return load_dict_contents_from_group(h5file, group)
        else:
            raise ValueError(group + ' not saved in ' + filename)

#end def load_dict_from_hdf5
#-----------------------------------------------


#-----------------------------------------------
def group_exist_in_hdf5(filename, group):
    """ Search if a group exist in a HDF5 file
    """

    group_exist = False
    with h5.File(filename, 'r') as h5file:
        if group in h5file.keys():
            group_exist = True

    return group_exist

#end def group_exist_in_hdf5
#-----------------------------------------------

