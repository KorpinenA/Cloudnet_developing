import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cloudnetpy.products.product_tools as ptools
from cloudnetpy import utils
import cloudnetpy.plotting as plot
from cloudnetpy.plot_meta import ATTRIBUTES

"""Plotting functions in development state
    - Generate figures to compare two variables
    - Generates figure of relative error of these two
"""


def _plot_relative_error(ax, error, axes, name):
    """Plots relative error"""
    pl = ax.pcolorfast(*axes, error[:-1, :-1].T, cmap='RdBu', vmin=-30,
                       vmax=30)
    colorbar = plot._init_colorbar(pl, ax)
    colorbar.set_label("%", fontsize=13)
    ax.set_title("Relative error of " + name, fontsize=14)


def generate_figure_from_two_files(nc_files, field_names, types, show=True,
                                   save_path=None, max_y=12, dpi=200):
    """ Compare product from two different files by subplotting
        Can plot several variables in loop,
        assuming all plotted variables being in one file

        field_names(list): list of plotted variables
        nc_files(list): list of two files to compare
                        [0] = old file
                        [1] = new file
    """
    # Not working if variable name not same for old and new file
    field_names = plot._parse_field_names(nc_files[0], field_names)
    data_fields = [ptools.read_nc_fields(nc_file, field_names) for nc_file in nc_files]

    # Information used for interpolation
    axes_data = [[plot._read_axes(nc_file, types[i]) for i in range(len(field_names))]
                  for nc_file in nc_files]
    data_fields, axes_data = _interpolate_select_variables(field_names, data_fields,
                                                           axes_data)
    case_date = plot._read_case_date(nc_files[0])
    subtit = (" from Cloudnet", " from CloudnetPy")

    for i, name in enumerate(field_names):
        fields = list(zip(*data_fields))[i]
        axes = list(zip(*axes_data))[i]

        plot_type = ATTRIBUTES[name].plot_type
        fig, ax = plot._initialize_figure(len(nc_files))
        for ii in range(len(nc_files)):
            field, axis = plot._fix_data_limitation(fields[ii], axes[ii], max_y)
            plot._set_axes(ax[ii], max_y)
            if plot_type == 'model':
                plot._plot_colormesh_data(ax[ii], field, name, axis)

            elif plot_type == 'bar':
                plot._plot_bar_data(ax[ii], field, name, axis)
                plot._set_axes(ax[ii], 1, ATTRIBUTES[name].ylabel)

            elif plot_type == 'segment':
                plot._plot_segment_data(ax[ii], field, name, axis)

            else:
                plot._plot_colormesh_data(ax[ii], field, name, axis)

            ax[ii].set_title(ATTRIBUTES[name].name + subtit[ii], fontsize=14)
        plot._add_subtitle(fig, len(nc_files), case_date)
        ax[-1].set_xlabel('Time (UTC)', fontsize=13)

        if save_path:
            file_name = parse_saving_name(save_path, case_date, [name], "_compare")
            plt.savefig(file_name, bbox_inches='tight', dpi=dpi)
        if show:
            plt.show()
        plt.close()


def generate_relative_err_fig(nc_files, field_names, types, show=True,
                              save_path=None, max_y=12, dpi=200):
    """ Creates figure of  relative error calculated two arrays with same size"""

    # Not working if variable name not same for old and new file
    field_names = plot._parse_field_names(nc_files[0], field_names)
    data_fields = [ptools.read_nc_fields(nc_file, field_names) for nc_file in nc_files]

    # Information used for interpolation
    axes_data = [[plot._read_axes(nc_file, types[i]) for i in range(len(field_names))]
                 for nc_file in nc_files]
    data_fields, axes_data = _interpolate_select_variables(field_names, data_fields,
                                                           axes_data)
    case_date = plot._read_case_date(nc_files[0])

    for i, name in enumerate(field_names):
        fields = list(zip(*data_fields))[i]
        axes = list(zip(*axes_data))[i]

        plot_type = ATTRIBUTES[name].plot_type
        fig, ax = plot._initialize_figure(1)

        plot._set_axes(ax[0], max_y)
        if plot_type == 'mesh':
            error = _calculate_relative_error(fields[0], fields[1])
            error, axis = plot._fix_data_limitation(error, axes[1], max_y)
            _plot_relative_error(ax[0], error, axis, name)

        if save_path:
            file_name = parse_saving_name(save_path, case_date, [name], "_error")
            plt.savefig(file_name, bbox_inches='tight', dpi=dpi)
        if show:
            plt.show()


def _interpolate_select_variables(field_names, data_field, axes):
    # Default is different dimensions with new and old files
    for i, name in enumerate(field_names):
        if ATTRIBUTES[name].plot_type == 'mesh':
            data_field[0][i] = _interpolate_data_and_dimensions\
                (data_field[0][i], list(zip(*axes))[i])
            axes[0][i] = axes[1][i]
    return data_field, axes


def parse_saving_name(save_path, case_date, field_names, ending):
    file_name = plot._create_save_name(save_path, case_date, field_names)
    file_name = file_name.split('.')
    file_name = file_name[0] + ending + ".png"
    return file_name


def _interpolate_data_and_dimensions(data, axes):
    n = np.min(data)
    data = np.asarray(data)
    data = utils.interpolate_2d(axes[0][0], axes[0][1], data,
                                axes[1][0], axes[1][1])
    data = ma.masked_where(data < n, data)
    return data


def _calculate_relative_error(old_data, new_data):
    ind = np.where((old_data > 0) & (new_data > 0))
    inds = np.full(new_data.shape, False, dtype=bool)
    inds[ind] = True
    old_data[~inds] = ma.masked
    new_data[~inds] = ma.masked
    return ((new_data - old_data) / old_data) * 100


if __name__ == '__main__':
    outname = '/home/korpinen/Documents/ACTRIS/cloudnetpy/test_data_iwc.nc'
    old_iwc_file = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_iwc-Z-T-method.nc'
    fname = '/home/korpinen/Documents/ACTRIS/cloudnet_data/categorize_test_file_new.nc'
    old_fname = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_categorize.nc'
    save_plots = '/home/korpinen/Documents/ACTRIS/cloudnetpy/plots/'
    """
    generate_figure_from_two_files([old_iwc_file, outname],
                                   ['iwc', 'iwc_error', 'iwc_retrieval_status'],
                                   [None, None, None])
    """
    generate_relative_err_fig([old_iwc_file, outname],
                              ['iwc', 'iwc_error', 'iwc_retrieval_status'],
                              [None, None, None])
