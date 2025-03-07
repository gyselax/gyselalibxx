# SPDX-License-Identifier: MIT

'''
General functions used for plot
'''

# pylint: disable=import-error

from io import BytesIO
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

import imageio


def plot_field1d(data_dict, titlename, filename, **opt_args):
    """
    Plots 1D field(s).

    Parameters
    ----------
    filename: pathlib.Path or str
        path of the file to generate

    data_dict : dictionary of the form
        { label : xarray.DataArray}
        the xarray is the 1D data to plot
        label to be put in the figure legend}
    """

    for data in data_dict.values():
        if data.ndim != 1:
            raise ValueError(f"Expected 1D data, got {data.ndim}D, {data.coords}")

    figsize = opt_args.get('figsize', (10, 10))
    fontsize = opt_args.get('fontsize', figsize[0]*2)
    fontname = opt_args.get('fontname', 'DejaVu Sans')
    title = opt_args.get('title', titlename)
    ymax = opt_args.get('ymax', None)
    ymin = opt_args.get('ymin', None)
    axis_font = opt_args.get(
        'axis_font', {'fontname': fontname, 'size': fontsize})
    title_font = opt_args.get('title_font', {'fontname': fontname,
                                             'size': fontsize,
                                             'verticalalignment': 'bottom'})
    scale = opt_args.get('scale', 'linear')
    linewidth = opt_args.get('linewidth', 2.)

    # Plotting with naive matplotlib
    mpl.style.use('classic')
    plt.rc('xtick', labelsize=fontsize*.8)
    plt.rc('ytick', labelsize=fontsize*.8)
    plt.rc('font', family=fontname)

    fig, ax = plt.subplots(figsize=figsize)
    axis_label = None
    for label, data in data_dict.items():
        ax.plot(data.coords[data.dims[0]], data.values, label=label, lw=linewidth)
        if axis_label is None:
            axis_label = data.dims[0]

    ax.set_xlabel(rf'${axis_label}$', **axis_font)
    ax.set_title(title, **title_font)
    ax.set_yscale(scale)
    ax.grid()
    ax.set_ylim(ymin, ymax)
    ax.margins(y=0.05)
    ax.axis('tight')
    ax.legend(loc='best')

    if isinstance(filename, str):
        filename = Path(filename)
    try:
        filename.parent.mkdir(parents=True)
    except FileExistsError:
        pass

    fig.savefig(filename)
    print('saved', filename)
    return fig, ax


def plot_field2d(data, titlename, filename, **opt_args):
    """
    Plots an image of a 2D field.

    Parameters
    ----------
    filename: pathlib.Path or str
        path of the file to generate

    data: xarray.DataArray
        2D data to plot
    """

    if data.ndim != 2:
        raise ValueError(f"Expected 2D data, got {data.ndim}D, {data.coords}")

    figsize = opt_args.get('figsize', (10, 10))
    fontsize = opt_args.get('fontsize', figsize[0]*2)
    fontname = opt_args.get('fontname', 'DejaVu Sans')
    title = opt_args.get('title', titlename)
    vmax = opt_args.get('vmax', np.nanmax(data))
    vmin = opt_args.get('vmin', np.nanmin(data))
    axis_font = opt_args.get(
        'axis_font', {'fontname': fontname, 'size': fontsize})
    title_font = opt_args.get('title_font', {'fontname': fontname,
                                             'size': fontsize,
                                             'verticalalignment': 'bottom'})
    scale = opt_args.get('scale', 'log')
    cmap = opt_args.get('cmap', 'jet')

    # Plotting with naive matplotlib
    mpl.style.use('classic')
    plt.rc('xtick', labelsize=fontsize*.8)
    plt.rc('ytick', labelsize=fontsize*.8)
    plt.rc('font', family=fontname)

    fig, ax = plt.subplots(figsize=figsize)
    if scale == 'log':
        if vmin == 0:
            vmin = 1e-10
            norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        elif vmin < 0 < vmax:
            print(f'negative values in {title}, with min={vmin}')
            vmin = 1e-10
            norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    im = ax.pcolormesh(data.coords[data.dims[1]],
                       data.coords[data.dims[0]],
                       data.values,
                       norm=norm,
                       shading='auto',
                       cmap=cmap)
    ax.set_title(title, **title_font)
    ax.set_xlabel(rf'${data.dims[1]}$', **axis_font)
    ax.set_ylabel(rf'${data.dims[0]}$', **axis_font)
    ax.axis('tight')
    fig.colorbar(im, ax=ax)

    if isinstance(filename, str):
        filename = Path(filename)
    try:
        filename.parent.mkdir(parents=True)
    except FileExistsError:
        pass

    fig.savefig(filename)
    print('saved', filename)
    return im, ax


def anim_field3d_append(data, out_writer, anim_dim='time', **opt_args):
    """
    Plots an animation of a 3D field over multiple images.

    The function can be used to add images to an already partially
    plotted animation.

    Parameters
    ----------
    out_writer: imageio.Writer
        image to append

    data: xarray.DataArray
        3D data to plot, the 1st

    anim_dim: str
        the dimension over which to generate the animation
    """
    if data.ndim != 3:
        raise ValueError("Expected 3D data, got {data.ndims}D, {data.coords}")

    for it in range(data.sizes[anim_dim]):
        it_val = float(data.coords[anim_dim][it])
        data_it = data.isel(indexers={anim_dim: it}, drop=True)
        with BytesIO() as imgdata:
            plot_field2d(data_it, filename=imgdata,
                         titlename=f'${anim_dim} = {it_val:02.2f}$',
                         **opt_args)
            imgdata.seek(0)
            out_writer.append_data(imageio.imread(imgdata))


def anim_field3d(data, out_file, **opt_args):
    """
    Plots an animation of a 3D field over multiple images.

    The function can be used to add images to an already partially
    plotted animation.

    Parameters
    ----------
    out_file: pathlib.Path or str
        path of the gif output file

    data: xarray.DataArray
        3D data to plot, the 1st

    anim_dim: str
        the dimension over which to generate the animation
    """
    if out_file.is_file():
        out_file.unlink()
    with imageio.get_writer(out_file, mode='I') as writer:
        anim_field3d_append(data, out_writer=writer, **opt_args)
