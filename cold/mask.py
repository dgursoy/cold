#!/usr/bin/env python3

import numpy as np
from scipy import signal, ndimage
import logging
import warnings
import xraydb
warnings.filterwarnings('ignore')


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'

from cold import core


PLANCK_CONSTANT = 6.58211928e-19  # [keV*s]
SPEED_OF_LIGHT = 299792458e+2  # [cm/s]


def wavelength(energy):
    """Return the wavelength [cm] for a given energy [keV]."""
    return 2 * np.pi * PLANCK_CONSTANT * SPEED_OF_LIGHT / energy


def discmask(geo, ind, inverted=True, exact=False, normalized=True, energy=10):
    # Mask create
    seq = np.load(geo['mask']['path'])
    if geo['mask']['reversed'] is True:
        seq = np.flip(seq)
    dseq = np.diff(seq)
    pt1 = np.zeros((64, 2))
    ind0 = 0
    ind1 = 0
    pointer = 0
    if seq[0] == 1:
        ind1 += 1
        pt1[0, 0] = pointer
    for m in range(seq.size - 1):
        pointer += geo['mask']['bitsizes'][seq[m]]
        if dseq[m] == 1: # we have a 1
            pt1[ind1, 0] = pointer
            ind1 += 1
        elif dseq[m] == -1: # we have a 0
            pt1[ind0, 1] = pointer
            ind0 += 1
    if seq[-1] == 1:
        pointer += geo['mask']['bitsizes'][seq[-1]]
        pt1[ind0, 1] = pointer
    if np.abs(geo['mask']['widening']) > 0:
        pt1[:, 0] -= geo['mask']['widening'] * 0.5
        pt1[:, 1] += geo['mask']['widening'] * 0.5

    # Rotation vector (intrinsic-yzx)
    alpha = geo['mask']['focus']['angley'] * np.pi / 180
    beta = geo['mask']['focus']['anglez'] * np.pi / 180
    gamma = geo['mask']['focus']['anglex'] * np.pi / 180
    rotmat = np.zeros((3, 3), dtype='float32')
    rotmat[0, 0] = np.cos(alpha) * np.cos(beta)
    rotmat[0, 1] = np.sin(alpha) * np.sin(gamma) - np.cos(alpha) * np.cos(gamma) * np.sin(beta)
    rotmat[0, 2] = np.cos(gamma) * np.sin(alpha) + np.cos(alpha) * np.sin(beta) * np.sin(gamma)
    rotmat[1, 0] = np.sin(beta)
    rotmat[1, 1] = np.cos(beta) * np.cos(gamma)
    rotmat[1, 2] = -np.cos(beta) * np.sin(gamma)
    rotmat[2, 0] = -np.cos(beta) * np.sin(alpha)
    rotmat[2, 1] = np.cos(alpha) * np.sin(gamma) + np.cos(gamma) * np.sin(alpha) * np.sin(beta)
    rotmat[2, 2] = np.cos(alpha) * np.cos(gamma) - np.sin(alpha) * np.sin(beta) * np.sin(gamma)

    # Rotation of mask axes
    mx = np.array([1, 0, 0], dtype='float32')
    my = np.array([0, 1, 0], dtype='float32')
    mz = np.array([0, 0, 1], dtype='float32')
    mx = np.dot(rotmat, mx)
    my = np.dot(rotmat, my)
    mz = np.dot(rotmat, mz)

    # Offset calculation
    p0 = core.pix2pos(ind, geo) # [<->, dis2det, v^]
    px = np.dot(p0, mx) 
    py = np.dot(p0, my) 
    pz = np.dot(p0, mz) 
    offset = np.abs(geo['mask']['thickness'] * px / py) 


    # Discretisize mask
    grid = creategrid(geo['mask'])
    mask = np.zeros(grid.size - 1)

    # Pad mask
    mask, grid = padmask(geo['mask'], mask, geo['mask']['pad'] / geo['mask']['resolution'])
    pt1[:, 0] += geo['mask']['pad']
    pt1[:, 1] += geo['mask']['pad']

    # Convolve
    for m in range(pt1.shape[0]):
        begin = int(np.ceil(pt1[m, 0] / geo['mask']['resolution']))
        end = int(np.floor(pt1[m, 1] / geo['mask']['resolution']))
        mask[begin+1:end+1] = 1
        mask[begin] = begin - pt1[m, 0] / geo['mask']['resolution']
        mask[end+1] = pt1[m, 1] / geo['mask']['resolution'] - end
    factor = 10
    mask = ndimage.zoom(mask, factor, order=1)
    if int(factor * offset / geo['mask']['resolution']) > 0:
        kernel = signal.tukey(int(factor * offset / geo['mask']['resolution']), alpha=0)
        kernel /= kernel.sum()
        mask = signal.convolve(mask, kernel, 'same')

    if exact == True:
        angpix = np.arctan(p0[0] / p0[1]) 
        angmsk = geo['mask']['focus']['anglez'] * np.pi / 180.
        mu = xraydb.mu_elam('Au', energy * 1e3) * 19.32 * geo['mask']['thickness'] * 1e-4 / np.cos(angpix + angmsk)
        mask = np.exp(-mu * mask)
        mask = core.invert(mask)

    mask = ndimage.zoom(mask, 1 / factor, order=1)

    if normalized == True:
        mask -= np.min(mask)
        mask /= np.max(mask)

    if inverted == True:
        mask = core.invert(mask)
    return mask



def mask(mask, inverted=True):
    """Returns a mask."""
    grid = creategrid(mask)
    vals = gridvals(mask, grid)
    vals, grid = padmask(mask, vals, mask['pad'] / mask['resolution'])
    if inverted == True:
        vals = core.invert(vals)
    return vals


def padmask(mask, vals, pad):
    vals = np.pad(vals, int(pad))
    padlen = 2 * pad * mask['resolution']
    totlen = masklength(mask) + padlen
    grid = np.arange(0, totlen, mask['resolution'])
    return vals, grid


def creategrid(mask):
    length = masklength(mask)
    grid = np.arange(0, length + mask['resolution'], mask['resolution'])
    return grid

def masklength(mask):
    sequence = loadmask(mask)
    n1 = numones(sequence)
    n0 = numzeros(sequence)
    length =  np.dot((n0, n1), mask['bitsizes'])
    return length


def plotmask(mask):
    """Plots the mask on a given grid."""
    import matplotlib.pyplot as plt
    plt.figure(figsize=(16, 1.5))
    plt.xlabel("Length [mu]")
    plt.ylabel("Mask")
    plt.plot(mask)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.close()


def loadmask(mask):
    return np.load(mask['path'])


def gridvals(mask, grid):
    sequence = np.load(mask['path'])
    if mask['reversed'] is True:
        sequence = np.flip(sequence)
    vals = np.zeros(grid.shape, dtype='float32')
    nbits = np.size(sequence)
    pointer = 0 
    for m in range(nbits):
        bit = sequence[m]
        size = mask['bitsizes'][sequence[m]]
        try: 
            start = grid[pointer]
            while (grid[pointer] - start < size):
                vals[pointer] = bit
                pointer += 1
        except IndexError:
            pass
    logging.info("Grid values assigned.")
    if mask['widening'] > 0:
        vals = widen(mask, vals)
    if mask['smoothness'] > 0:
        vals = smooth(mask, vals, mask['alpha'])
    return vals


def widen(mask, vals):
    size = int(mask['widening'] / mask['resolution'])
    kernel = signal.tukey(size, alpha=0.0)
    kernel /= kernel.sum()
    vals = signal.convolve(vals, kernel, 'same')
    vals[vals > 1e-6] = 1
    logging.info("Mask widened.")
    return vals


def smooth(mask, vals, alpha):
    size = int(mask['smoothness'] / mask['resolution'])
    kernel = signal.tukey(size, alpha=alpha)
    kernel /= kernel.sum()
    vals = signal.convolve(vals, kernel, 'same')
    logging.info("Mask smoothed.")
    return vals


def numbits(sequence):
    return np.size(sequence)


def numones(sequence):
    return np.sum(sequence)


def numzeros(sequence):
    return numbits(sequence) - numones(sequence)

def diffvals(mask):
    return -np.diff(mask)
