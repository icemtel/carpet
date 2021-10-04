'''
Functions to pickle/unpickle

Functions to load/save data from time-series:
- one trajectory may be split into parts

'''
import os
import numpy as np


def define_dump_object(default_output_folder='.'):
    '''
    Define a function-shortcut for pickling objects.
    - Pickling is a way to save python objects to the disk space.
    - complex objects might get broken (e.g. matplotlib figure)
    - works well with dictionaries, numpy.arrays
    - There is a chance that won't be able to load objects on a different python/libraries version.
    '''
    import pickle, os

    def dump_object(obj, filename, path=default_output_folder):
        filename = os.path.join(path, filename)
        print(filename)
        with open(filename, 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

    return dump_object


def define_load_object(default_output_folder='.'):
    '''
    Define a function-shortcut for loading (un-pickling) objects.
    - There is a chance that won't be able to load objects on a different python/libraries version.
    - defines "load_object" function with the given default folder
    '''
    import pickle, os

    def load_object(filename, path=default_output_folder):
        filename = os.path.join(path, filename)
        with open(filename, 'rb') as f:
            obj = pickle.load(f)
        return obj

    return load_object


def define_get_filename(base_filename):
    def get_filename(irun, ipart, path):
        if ipart == 0:
            return os.path.join(path, base_filename + f'_{irun}.npy')
        else:
            return os.path.join(path, base_filename + f'_{irun}_pt{ipart}.npy')

    return get_filename


def define_loader(base_filename, npart_max=64):
    get_filename = define_get_filename(base_filename)

    def load_xs(irun, path='.', ncycle=None):
        '''
        Read saved file with xs array
        :param ncycle: if given cut trajectory up to a certain number of cycles;
        '''
        xs_list = []

        for ipart in range(npart_max):
            filename = get_filename(irun, ipart, path)
            if os.path.isfile(filename):
                xs_part = np.load(filename)
                xs_list.append(xs_part)
            else:
                break
        xs = np.concatenate(xs_list)  # trajectory glued back from parts
        if ncycle is not None:
            xs = xs[:ncycle]

        return xs

    return load_xs


def define_load_ts(base_filename='ts', npart_max=64):
    '''
    With timestamps it's a bit different - in every part of the trajectory my times start from 0.
    With this function I take this into account
    :param name:
    :return:
    '''
    get_ts_filename = define_get_filename(base_filename)

    def load_ts(irun, path='.', ncycle=None):
        '''
        Read pickled file with phis array
        :param ncycle: if given cut trajectory up to a certain number of cycles;
        '''
        ts_list = []

        for ipart in range(npart_max):
            filename = get_ts_filename(irun, ipart, path)
            if os.path.isfile(filename):
                ts_part = np.load(filename)
                if ipart > 0:
                    ts_part += t_last  # add up previous time (in each part time starts from zero)
                ts_list.append(ts_part)
                t_last = ts_part[-1]
            else:
                break
        ts = np.concatenate(ts_list)  # trajectory glued back from parts
        # try:
        #     phis = np.concatenate(phis_list)  # trajectory glued back from parts
        # except Exception as e:
        #     print('FAIL; params:', D, dt, irun, ncycle)
        #     raise e
        if ncycle is not None:
            ts = ts[:ncycle]
        return ts

    return load_ts

if __name__ == '__main__':
    load_phis = define_loader('phi')
    load_dphidt_mean = define_loader('dphidt_mean')
    load_dphidt_std = define_loader('dphidt_std')
    get_ts_filename = define_get_filename('ts')