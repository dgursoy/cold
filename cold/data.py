#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import re
import os
import yaml
import logging
import warnings
warnings.filterwarnings('ignore')


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'



def config(path):
    stream = open(path, 'r')
    pars = yaml.safe_load(stream)
    return pars['file'], pars['geo'], pars['algo']



def loadsingle(file, id=0):
    """Loads Laue diffraction data."""
    if file['stacked'] is True:
        files = loadstack(file)
        if file['ext'] == 'h5':
            values = loadh5(files[id], file['h5']['key'])  
    else:
        if file['ext'] == 'h5':
            values = loadh5(file['path'], file['h5']['key'])[id:id+1]
    return values


def load(file, collapsed=True, partitioned=True, index=None):
    """Loads Laue diffraction data."""
    if file['stacked'] is True:
        files = loadstack(file)
        if file['ext'] == 'h5':
            values = loadh5files(files, file['h5']['key'])
    else:
        if file['ext'] == 'h5':
            begin, end = file['range']
            values = loadh5(file['path'], file['h5']['key'])[begin:end]
            values = np.swapaxes(values, 0, 2)
            values = np.swapaxes(values, 0, 1)
            values = values.copy()
    if index is None:
        index = cherrypickpixels(values, file['threshold'], file['frame'])
    if collapsed is True:
        values = collapse(values, index)
    if partitioned is True:
        values = partition(values, file['chunks'])
        index = partition(index, file['chunks'])
    return values, index


def partition(data, n):
    k, m = divmod(data.shape[0], n)
    return list(data[i * k + min(i, m):
        (i + 1) * k + min(i + 1, m)] for i in range(n))


def loadstack(file):
    """Loads Laue diffraction data stack."""
    files = getfiles(file['path'], file['ext'])
    files = cherrypickfiles(files, file['range'])
    return files


def getfiles(path, ext):
    """Returns the sorted list of file names with a given 
    extension in a directory."""
    files = Path(path).glob('*.' + ext)
    return sorted(files, 
        key=lambda path: int(re.sub('\D', '', str(path))))


def cherrypickfiles(files, range=None):
    begin, end = range
    if range is None:
        range = [0, len(files)]
    return files[begin:end]


def loadh5files(files, key):
    data = initdata(files, key)
    for m in range(len(files)):
        data[:, :, m] = loadh5(files[m], key)
    return data


def initdata(files, key):
    img = loadh5(files[0], key)
    nx, ny = img.shape
    return np.zeros((nx, ny, len(files)), dtype='float32')


def loadh5(file, key):
    import h5py
    f = h5py.File(file, 'r')
    value = f[key][:]
    logging.info("Loaded: " + str(file))
    return value


def save(path, data, swap=False):
    """Saves Laue diffraction data to file."""
    import dxchange
    if swap is True:
        data = np.swapaxes(data, 0, 2)
        data = np.swapaxes(data, 1, 2)
    dxchange.write_tiff(data, fname=path)
    logging.info("Saved: " + str(path) + ".tiff")


def cherrypickpixels(data, threshold, frame):
    """Returns pixel indices above a threshold and inside a frame."""
    pixels = tagpixels(data, threshold)
    pixels = rejectpixels(frame, pixels)
    index = getindices(pixels)
    return index


def tagpixels(data, threshold):
    """Return a binary image that shows tagged pixels."""
    img = np.mean(data, axis=2)
    pixels = np.zeros(img.shape, dtype='int16')
    pixels[img >= threshold] = 1
    pixels[img == 65535] = 0
    return pixels


def rejectpixels(frame, pixels):
    """Return an image by zeroing pixels outside a frame."""
    pixels[0:frame[0], :] = 0
    pixels[frame[1]::, :] = 0
    pixels[:, 0:frame[2]] = 0
    pixels[:, frame[3]::] = 0
    return pixels


def getindices(pixels):
    """Computes the indices for tagged pixels from a binary image."""
    ix, iy = np.where(pixels == 1)
    index = np.array([ix, iy], dtype='int32').T
    return index


def collapse(data, index):
    """Collapses an ND array into a (N-1)D array."""
    numpixels = index.shape[0]
    shape = np.array(data.shape)
    shape[1] = numpixels
    _data = np.zeros(shape[1:], dtype='float32')
    for m in range(numpixels):
        _data[m] = data[index[m, 0], index[m, 1]]
    return _data


def expand(data, index, shape):
    """Expands an array."""
    nx = index.shape[0]
    try:
        dx = data.shape[1]
    except IndexError:
        dx = 1    
    _data = np.zeros((shape[0], shape[1], dx), dtype='float32')
    for m in range(nx):
        _data[index[m, 0], index[m, 1]] = data[m]
    return _data.squeeze()


def pack(data):

    # Number pf processes
    chunks = len(data)

    # Number of total pixels
    npix = 0
    for m in range(chunks):
        npix += data[m].shape[0]
    
    # Initialize array
    try:
        dy = data[m].shape[1]
        arr = np.zeros((npix, dy), dtype=data[0].dtype)
    except IndexError:
        arr = np.zeros((npix, ), dtype=data[0].dtype) 

    # Write to array
    ind = 0
    for m in range(chunks):
        begin = ind
        end = begin + data[m].shape[0]
        ind += data[m].shape[0]
        arr[begin:end] = data[m]
    return arr


def saveimg(path, vals, inds, shape, swap=False):
    _vals = packandexpand(vals, inds, shape)
    save(path, _vals, swap)


def packandexpand(vals, inds, shape):
    _vals = pack(vals)
    _inds = pack(inds)
    _vals = expand(_vals, _inds, shape)
    return _vals


def saveplt(path, vals, grid):
    import matplotlib.pyplot as plt
    p = Path(path).parents[0]
    if not os.path.exists(p):
        os.makedirs(p)
    _vals = sum(vals)
    _grid = np.arange(*grid)
    plt.figure(figsize=[20, 5])
    plt.subplot(211)
    plt.step(_grid, _vals)
    plt.grid()
    plt.subplot(212)
    _vals[_vals < 1] = 0
    plt.semilogy(_grid, _vals, drawstyle='steps')
    plt.grid()
    plt.xlabel('[mm]')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    logging.info("Saved: " + str(path) + ".png")


def plotarr(path, arr, plots=False):
    import matplotlib.pyplot as plt
    import dxchange
    arr = pack(arr)
    dxchange.write_tiff(arr, path)
    logging.info("Saved: " + str(path) + ".tiff")
    if plots is True:
        p = Path(path).parents[0]
        if not os.path.exists(p):
            os.makedirs(p)
        for m in range(len(arr)):
            plt.figure(figsize=(8, 3))
            plt.plot(arr[m], 'ro-')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(str(path) + "-" + str(m))
            plt.close()
            logging.info("Saved: " + str(path) + "-" + str(m) + ".png")