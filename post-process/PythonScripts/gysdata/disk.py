# SPDX-License-Identifier: MIT

'''
Functions used to create a store that respresents the outputs
from the code on disk
'''

# pylint: disable=import-error

from collections import OrderedDict
from pathlib import Path
import re
from yaml.loader import SafeLoader
import yaml

import dask.array as da
import h5py as h5
import xarray as xr

from . import coord_getters


class FieldInfo(object):
    '''
    Get field info
    '''
    def __init__(self, dimensions=None, path=None):
        # an OrderedDict of coord_getters
        self.dimensions = dimensions if dimensions is not None else OrderedDict()
        # a list of DsetInfo
        self.path = path if path is not None else []

    def __str__(self):
        return ("( path: [" + ", ".join([str(p) for p in self.path]) + "], dimensions:{"
                + ", ".join([k+": "+str(v) for k, v in self.dimensions.items()])+"})")


class DsetInfo(object):
    '''
    Set info
    '''
    def __init__(self, file=None, dset=None):
        self.file = file
        self.dset = dset

    def __str__(self):
        return "( file: "+str(self.file)+" dset: "+str(self.dset)+")"


class Store(object):
    '''
    A Store that respresents the outputs from the code on disk.
    '''

    def __init__(self, data_dir, data_structure=None, fields=None):
        '''
        Build a new Store

        Parameters
        ----------
        data_dir: pathlib.Path or str
            Directory that holds the data described in the data_structure file.

        data_structure: pathlib.Path or str
            Yaml file that describes the files content.

        fields: str or list of str or None, optional
            List of all the fields that will be built.
            By default, all the fields described in data_structure are built.

        Returns
        -------
        self: Store instance
        '''

        # dictionary of name to DataArray
        self.__datadict = {}

        # The base data dir
        data_dir = Path(data_dir).resolve()

        if data_structure is None:
            data_structure = Path('data_structure.yaml')
        if not data_structure.is_absolute():
            if (data_dir / data_structure).is_file():
                data_structure = data_dir/data_structure
            else:
                src_dir = Path(__file__).absolute().parent
                data_structure = src_dir/data_structure

        # Fields structure as a dictionary: dict[str, FieldInfo]: field_name -> description
        fields_structure = {}
        if fields is not None and not isinstance(fields, list):
            self.fields = [fields]
        with open(data_structure) as f:
            fields_structure_data = yaml.load(f, Loader=SafeLoader)
            for field_name, field_data in fields_structure_data.items():
                if fields is not None and field_name not in fields:
                    continue

                field_infos = fields_structure.setdefault(
                    field_name, FieldInfo())

                for dimension_data in field_data['dimensions']:
                    for dimension_name, dimension_info in dimension_data.items():
                        coord_getter_name, coord_getter_args = list(
                            dimension_info.items())[0]
                        coord_getter_class = getattr(
                            coord_getters, coord_getter_name)
                        coord_getter = coord_getter_class(data_dir,
                                                          *coord_getter_args)
                        field_infos.dimensions[dimension_name] = coord_getter

                path_list = field_data['path']
                if not isinstance(path_list, list):
                    path_list = [path_list]
                for path in path_list:
                    file_re = re.compile(path['file'])
                    dset_re = re.compile(path['dataset'])
                    field_infos.path.append(DsetInfo(file_re, dset_re))

        # Fields in each file: dict[str, list[str]]: file_name -> fields_names
        file_to_fields = {}
        for file in data_dir.glob('**/*'):
            if file.is_file():
                for field_name, field_info in fields_structure.items():
                    for path in field_info.path:
                        relfile = file.relative_to(data_dir)
                        if path.file.fullmatch(str(relfile)) is not None:
                            fields_list = file_to_fields.setdefault(
                                file.resolve(), [])
                            fields_list.append(field_name)

        # List that keeps track of which coordinate have already been read to prevent reading the
        # same datasets multiples times
        global_coords = {}

        # Iterate over all files first to open them as few times as possible
        for path, ffields in file_to_fields.items():
            h5file = h5.File(path, 'r')
            for field_name in ffields:
                field = fields_structure[field_name]
                for dset_info in field.path:
                    file_match = dset_info.file.fullmatch(
                        str(path.relative_to(data_dir)))
                    if file_match is None:
                        continue
                    for dset_name in h5file.keys():
                        dset_match = dset_info.dset.fullmatch(dset_name)
                        if dset_match is None:
                            continue
                        coords = OrderedDict()
                        for dim_name, coord_getter in field.dimensions.items():
                            coords[dim_name] = coord_getter(
                                field,
                                file_match,
                                h5file,
                                [],
                                dset_name,
                                dset_match,
                                global_coords
                            )

                        data = da.from_array(h5file[dset_name])

                        # Small verification the match between the coordinates and the data.
                        target_ndim = len(coords)
                        if data.ndim > target_ndim:
                            raise ValueError("Not enough dimension specified for dataset "
                                             f"{dset_name} in {path}: expected {data.ndim} "
                                             f"dimensions, got {coords}.")
                        if data.ndim < target_ndim:
                            new_shape = (1, ) * (target_ndim -
                                                 data.ndim) + data.shape
                            data = data.reshape(new_shape)

                        kv_coords = list(coords.items())
                        for dim in range(data.ndim):
                            if data.shape[dim] < len(kv_coords[dim][1]):
                                coords[kv_coords[dim][0]] = kv_coords[dim][1][0:data.shape[dim]]
                        try:
                            data_array = xr.DataArray(
                                data, dims=coords.keys(), coords=coords)
                        except Exception as e:
                            raise ValueError(
                                str(e)+" for '{}'".format(field_name))
                        # remove the names to work-around xr.combine_by_coords failure,
                        # see https://github.com/pydata/xarray/pull/5834
                        data_array.name = None

                        self.__datadict.setdefault(
                            field_name, []).append(data_array)

        # Regroup DataArrays by coords to have one DataArray per field.
        for k, v in self.__datadict.items():
            self.__datadict[k] = xr.combine_by_coords(v)
            self.__datadict[k].name = k

    def get_data(self):
        '''
        Return the dictionary associated to the data structure
        '''
        return self.__datadict

    def __getitem__(self, field):
        '''
        return an xarray.DataArray (or a dataset if field is a list)
        '''
        return self.__datadict[field]
