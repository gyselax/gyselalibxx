# SPDX-License-Identifier: MIT

'''
Functions used to define the coordinates
'''

import re

import h5py as h5
import xarray as xr

class local_coord:
    '''
    Local coordinates
    '''

    def __init__(self, data_dir, key):
        self.key = key

    def __call__(self,
                 _,
                 __,
                 file,
                 ___,
                 ____,
                 _____,
                 ______,):
        return [file[self.key][()]]


class constant_coord:
    '''
    Constant coordinates
    '''

    def __init__(self, data_dir, *values):
        self.values = list(values)

    def __call__(self,
                 _,
                 __,
                 ___,
                 ____,
                 _____,
                 ______,
                 _______):
        return self.values


class global_coord:
    '''
    Global coordinates
    '''

    def __init__(self, data_dir, filename, key):
        self.filename = data_dir / filename
        self.key = key
    def __call__(self,
                 _,
                 __,
                 ___,
                 ____,
                 _____,
                 ______,
                 global_coords):
        try:
            coords = global_coords[self.filename / self.key]
            return coords
        except KeyError:
            file = h5.File(self.filename, 'r')
            tmp = file[self.key][...]
            global_coords[self.filename / self.key] = tmp
            file.close()
            return tmp



class filename_coord:
    '''
    Local coordinates
    '''

    def __init__(self, data_dir, key):

        self.key = key
    def __call__(self,
                 _,
                 regexp_match,
                 __,
                 ___,
                 ____,
                 _____,
                 ______):
        return [int(regexp_match.group(self.key))]



class scattered_coord:
    '''
    Scattered coordinates
    '''

    def __init__(self,
                 data_dir,
                 scatter_pattern_file_regexp,
                 scatter_pattern_dataset_regexp,
                 scatter_pattern_index,
                 filename_regexp_group_key,
                 global_coords_key,
                 global_coords_file,
                 ):

        # Where to find the pattern in 3 steps :
        # 1. which file ?
        self.scatter_pattern_file_regexp = re.compile(scatter_pattern_file_regexp)
        # 2. which dataset in said file ?
        self.scatter_pattern_dataset_regexp = re.compile(scatter_pattern_dataset_regexp)
        # 3. where in that dataset ?
        self.scatter_pattern_index = scatter_pattern_index

        # Now that we have the pattern where are we located in it ?
        self.filename_regexp_group_key = filename_regexp_group_key

        # Ok, and which dimension are we looking for ?
        self.global_coords_key = global_coords_key
        # If it hasn't already been logged, where do we find it ?
        self.global_coords_file = data_dir + global_coords_file



    def __call__(self,
                 field,
                 file_match,
                 _,
                 filelist,
                 __,
                 dataset_match,
                 global_coords):
        try:
            pattern = global_coords[field+"pattern"]
        except KeyError:
            pattern = self.build_pattern(filelist)
            global_coords[field+"pattern"] = pattern

        file_match_dictionnary = file_match.groupdict()
        dataset_match_dictionnary = dataset_match.groupdict()

        match_dictionnary = {}

        match_dictionnary.update(file_match_dictionnary)
        match_dictionnary.update(dataset_match_dictionnary)

        for dim in pattern.dims:
            if dim in match_dictionnary.keys():
                pattern = pattern[{dim: int(match_dictionnary[dim])}]
        pattern = pattern.data

        # This is the number of fragments in which the dimension is scattered
        n_fragments_dim = pattern[self.scatter_pattern_index]

        # This is the fragment number we are looking for
        fragment_index = int(match_dictionnary[self.filename_regexp_group_key])

        try:
            coords_full = global_coords[self.global_coords_file + self.global_coords_key]
        except KeyError:
            coords_file = h5.File(self.global_coords_file, 'r')
            coords_full = coords_file[self.global_coords_key][...]
            global_coords[self.global_coords_key] = coords_full
            coords_file.close()

        if len(coords_full)% 2:
            coords_full = coords_full[:-1]
        i_start = fragment_index * (len(coords_full) // n_fragments_dim)
        i_end = i_start + (len(coords_full) // n_fragments_dim)

        return coords_full[i_start:i_end]

    def build_pattern(self, filelist):
        '''
        Build the required pattern
        '''
        dataArray_list = []

        relevant_files = []
        # looking through the filelist to find the relevant ones
        for f in filelist:
            match = self.scatter_pattern_file_regexp.search(f)
            if match is not None:
                relevant_files.append((f, match.groupdict()))

        for f, file_groups in relevant_files:
            file = h5.File(f, 'r')
            relevant_dsets_keys = []

            #finding the relevant dataset_keys
            for k in file.keys():
                match_key = self.scatter_pattern_dataset_regexp.search(k)
                if match_key is not None:
                    relevant_dsets_keys.append((k, match_key.groupdict()))

            #build the dataArray that holds the pattern.
            for key, key_groups in relevant_dsets_keys:
                data = file[key][...]
                dims = ["dim_"+str(i) for i in range(data.ndim)]
                coords_dict = {}

                #Labeling dimensions
                if len(file_groups) != 0:
                    dims = dims + list(file_groups.keys())
                    coords_dict.update(file_groups)
                    data = data.reshape(data.shape + (1,) * len(file_groups))

                if len(key_groups) != 0:
                    dims = dims + list(key_groups.keys())
                    coords_dict.update(key_groups)
                    data = data.reshape(data.shape + (1,) * len(key_groups))

                for k, v in coords_dict.items():
                    if not isinstance(v, list):
                        coords_dict[k] = [int(v)]

                dataArray = xr.DataArray(data, dims=dims, coords=coords_dict)
                dataArray_list.append(dataArray)
            file.close()
        pattern = xr.combine_by_coords(dataArray_list)

        return pattern
