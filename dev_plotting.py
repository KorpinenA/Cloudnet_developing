import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cloudnetpy.products.product_tools as ptools
from cloudnetpy import utils
import cloudnetpy.plotting as plot
from cloudnetpy.plot_meta import ATTRIBUTES
from cloudnetpy.meta_for_old_files import fix_old_data

"""Plotting functions in development state
    - Generate figures to compare two variables
"""


def generate_figure_from_two_files(nc_files, field_names, types, old_file=False, show=True,
                                   re_err=False, save_path=None, max_y=12, dpi=200):
    """ Compare product from two different files by subplotting
        Can plot several variables in loop,
        assuming all plotted variables being in one file

        field_names(list): list of plotted variables
        nc_files(list): list of two files to compare
                        [0] = old file
                        [1] = new file
        types(list): list of plot types of variables, None if not Model
    """
    old_data, new_data = [plot._find_valid_fields(nc_file, field_names) for nc_file in nc_files]
    data_fields = [old_data[0], new_data[0]]
    field_names = new_data[-1]

    # Information used for interpolation
    axes_data = [[plot._read_axes(nc_file, types[i]) for i in range(len(field_names))]
                  for nc_file in nc_files]
    data_fields, axes_data = _interpolate_select_variables(field_names, data_fields,
                                                           axes_data)
    case_date, site_name = plot._read_case_date(nc_files[0])
    subtit = (" from Cloudnet", " from CloudnetPy")

    for i, name in enumerate(field_names):
        fields = list(zip(*data_fields))[i]
        axes = list(zip(*axes_data))[i]

        plot_type = ATTRIBUTES[name].plot_type
        if plot_type == 'mesh' and re_err:
            fig, ax = plot._initialize_figure(len(nc_files) + 1)
        else:
            fig, ax = plot._initialize_figure(len(nc_files))

        for ii in range(len(nc_files)):
            field, axis = plot._fix_data_limitation(fields[ii], axes[ii], max_y)
            plot._set_axes(ax[ii], max_y)
            if plot_type == 'model':
                plot._plot_colormesh_data(ax[ii], field, name, axis)

            elif plot_type == 'bar':
                if name == 'lwp' and ii == 0:
                    field = field*1000
                plot._plot_bar_data(ax[ii], field, name, axis[0])
                plot._set_axes(ax[ii], 1, ATTRIBUTES[name].ylabel)

            elif plot_type == 'segment':
                if old_file and ii == 0:
                    field, name = fix_old_data(field, name)
                plot._plot_segment_data(ax[ii], field, name, axis)

            else:
                plot._plot_colormesh_data(ax[ii], field, name, axis)
                if re_err and ii == len(nc_files) - 1:
                    plot._set_axes(ax[-1], max_y)
                    error = _calculate_relative_error(fields[0], fields[1])
                    error, axis = plot._fix_data_limitation(error, axes[1], max_y)
                    _plot_relative_error(ax[-1], error, axis, name)

            ax[ii].set_title(ATTRIBUTES[name].name + subtit[ii], fontsize=14)
        plot._add_subtitle(fig, case_date, site_name)
        ax[-1].set_xlabel('Time (UTC)', fontsize=13)

        if save_path:
            file_name = _parse_saving_name(save_path, case_date, max_y, [name], "_compare")
            plt.savefig(file_name, bbox_inches='tight', dpi=dpi)
        if show:
            plt.show()
        plt.close()


def _interpolate_select_variables(field_names, data_field, axes):
    # Default is different dimensions with new and old files
    def _interpolate_data_and_dimensions(data, axes):
        n = np.min(data)
        data = np.asarray(data)
        data = utils.interpolate_2d(axes[0][0], axes[0][1], data,
                                    axes[1][0], axes[1][1])
        data = ma.masked_where(data < n, data)
        return data
    
    for i, name in enumerate(field_names):
        if ATTRIBUTES[name].plot_type == 'mesh':
            data_field[0][i] = _interpolate_data_and_dimensions\
                (data_field[0][i], list(zip(*axes))[i])
            axes[0][i] = axes[1][i]
    return data_field, axes


def _calculate_relative_error(old_data, new_data):
    ind = np.where((old_data >= np.min(old_data)) & (new_data >= np.min(new_data)))
    inds = np.full(new_data.shape, False, dtype=bool)
    inds[ind] = True
    old_data[~inds] = ma.masked
    new_data[~inds] = ma.masked
    return ((new_data - old_data) / old_data) * 100


def _plot_relative_error(ax, error, axes, name):
    """Plots relative error"""
    pl = ax.pcolorfast(*axes, error[:-1, :-1].T, cmap='RdBu', vmin=-30,
                       vmax=30)
    colorbar = plot._init_colorbar(pl, ax)
    colorbar.set_label("%", fontsize=13)
    ax.set_title("Relative error of " + name, fontsize=14)


def _parse_saving_name(save_path, case_date, max_y, field_names, ending):
    file_name = plot._create_save_name(save_path, case_date, max_y, field_names)
    file_name = file_name.split('.')
    file_name = file_name[0] + ending + ".png"
    return file_name


def storage_file_info(quantity):
    cloudnetpy_class = '/home/korpinen/Documents/ACTRIS/cloudnetpy/test_data_classification.nc'
    cloudnetpy_iwc = '/home/korpinen/Documents/ACTRIS/cloudnetpy/test_data_iwc.nc'
    cloudnetpy_lwc = '/home/korpinen/Documents/ACTRIS/cloudnetpy/test_data_lwc.nc'
    cloudnetpy_drizzle = '/home/korpinen/Documents/ACTRIS/cloudnetpy/test_data_drizzle.nc'

    cloudnet_class = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_classification.nc'
    cloudnet_iwc = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_iwc-Z-T-method.nc'
    cloudnet_lwc = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_lwc-scaled-adiabatic.nc'
    cloudnet_drizzle = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_drizzle.nc'

    cloudnetpy_cat = '/home/korpinen/Documents/ACTRIS/cloudnet_data/categorize_test_file_new.nc'
    cloudnet_cat = '/home/korpinen/Documents/ACTRIS/cloudnet_data/20181204_mace-head_categorize.nc'

    class_quantities = ['target_classification', 'detection_status']
    iwc_quantities = ['iwc', 'iwc_inc_rain', 'iwc_error', 'iwc_retrieval_status']
    lwc_quantities = ['lwc', 'lwc_error', 'lwc_retrieval_status', 'lwp']
    drizzle_quantities = ['Do', 'mu', 'S', 'N', 'drizzle_N', 'lwc', 'drizzle_lwc',
                          'lwf', 'drizzle_lwf', 'drizzle_retrieval_status',
                          'droplet_fall_velocity', 'v_drizzle', 'w', 'v_air']
    drizzle_errors = ['Do_error', 'S_error', 'N_error', 'drizzle_N_error',
                      'lwc_error', 'drizzle_lwc_error', 'lwf_error', 'drizzle_lwf_error',
                      'droplet_fall_velocity_error', 'v_drizzle_error']
    cat_quantities = ['Z', 'Z_error', 'beta', 'Ze', 'ldr', 'width', 'v', 'lwp',
                      'radar_liquid_atten', 'radar_gas_atten']
    cat_bits = ['droplet', 'falling', 'cold', 'melting', 'aerosol', 'insect']
    qual_bits = ['radar', 'lidar', 'clutter', 'molecular', 'attenuated', 'corrected']

    data = {'classification': [cloudnet_class, cloudnetpy_class, class_quantities, [None] * len(class_quantities)],
            'iwc': [cloudnet_iwc, cloudnetpy_iwc, iwc_quantities, [None] * len(iwc_quantities)],
            'lwc': [cloudnet_lwc, cloudnetpy_lwc, lwc_quantities, [None] * len(lwc_quantities)],
            'drizzle': [cloudnet_drizzle, cloudnetpy_drizzle, drizzle_quantities, [None] * len(drizzle_quantities)],
            'drizzle_error': [cloudnet_drizzle, cloudnetpy_drizzle, drizzle_errors, [None] * len(drizzle_errors)],
            'categorize': [cloudnet_cat, cloudnetpy_cat, cat_quantities, [None] * len(cat_quantities)],
            'cat_bits': [cloudnet_cat, cloudnetpy_cat, cat_bits, [None] * len(cat_bits)],
            'qual_bits': [cloudnet_cat, cloudnetpy_cat, qual_bits, [None] * len(qual_bits)]}

    return data[quantity]


if __name__ == '__main__':
    # Run program of compering old and new files
    save_plots = '/home/korpinen/Documents/ACTRIS/cloudnetpy/plots/'

    # Key options: 'classification','iwc', 'lwc', 'drizzle', 'drizzle_error', 'categorize', 'cat_bits', 'qual_bits'
    data = storage_file_info('lwc')
    generate_figure_from_two_files([data[0], data[1]],
                                   data[2], data[3],
                                   old_file=True, max_y=12)
