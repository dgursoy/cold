#!/usr/bin/env python3

import numpy as np
import ctypes
import cv2
import logging
from cold import recpos, recsig, rectau, recdec, footprint, convolve, pack, saveplt
import multiprocessing
import os


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'


LIBCOLD = ctypes.cdll.LoadLibrary(
    'build/lib.macosx-10.7-x86_64-3.6/libcold.cpython-36m-darwin.so')



def decode(data, mask, grid, geo, algo):
    """Returns the mask position and the signal."""
    data = data.copy()
    mask = np.abs(mask.copy() - 1)
    ratio = getratio(grid, geo['mask']['step'])
    for m in range(len(data)):
        data[m] = stretch(data[m], ratio)
    msk = [mask] * len(data)
    pos = [np.array(algo['init'][0] * ratio, dtype='int32')] * len(data)
    sig = [np.array(algo['init'][1] * ratio, dtype='int32')] * len(data)
    tau = [np.array(algo['init'][2], dtype='float32')] * len(data)
    dec = [np.array(algo['init'][3] / np.sqrt(ratio), dtype='float32')] * len(data)
    alg = [algo] * len(data)
    rat = [ratio] * len(data)
    results = runpar(_decode, [data, msk, pos, sig, tau, dec, alg, rat], len(data))
    for m in range(len(data)):
        pos[m] = results[m][0] / ratio
        sig[m] = results[m][1] / ratio
        tau[m] = results[m][2]
        dec[m] = results[m][3] * np.sqrt(ratio)
    return pos, sig, tau, dec


def _decode(args):
    """Returns the mask position and the signal."""
    data = args[0]
    mask = args[1]
    pos0 = args[2]
    sig0 = args[3]
    tau0 = args[4]
    dec0 = args[5]
    algo = args[6]
    ratio = args[7]
    npix = data.shape[0]
    pos = np.zeros((npix, ), dtype='int32')
    sig = np.zeros((npix, ), dtype='int32')
    tau = np.zeros((npix, ), dtype='float32')
    dec = np.zeros((npix, ), dtype='float32')
    for m in range(data.shape[0]):
        pos[m], sig[m], tau[m], dec[m] = pixsolve(data[m], mask, pos0, sig0, tau0, dec0, algo)
        if pos[m] > 0:
            logging.info("Pixel decoded: " + 
                str(m) + '/' + str(npix - 1) + 
                " pos=" + str(pos[m] / ratio) + 
                " sig=" + str(sig[m] / ratio) + 
                " tau=" + str(tau[m]) + 
                " dec=" + str(dec[m] * np.sqrt(ratio)))
    return pos, sig, tau, dec


def pixsolve(data, mask, pos, sig, tau, dec, algo):
    """Returns the mask position and the signal footprintof  a single pixel."""
    for m in range(algo['iters']):
        pos = recpos(data, mask, pos, sig, tau, dec, algo['method'])
        sig = recsig(data, mask, pos, sig, tau, dec, algo['method'])
        tau = rectau(data, mask, pos, sig, tau, dec, algo['method'])
        dec = recdec(data, mask, pos, sig, tau, dec, algo['method'])
    return pos, sig, tau, dec


def stackargs(args, nproc):
    argstack = []
    for m in range(nproc):
        inputs = []
        for arg in args:
            inputs.append(arg[m])
        argstack.append(inputs)
    return argstack


def runpar(myfunc, args, nproc):
    pool = multiprocessing.Pool(processes=nproc)
    results = pool.map(myfunc, stackargs(args, nproc))
    pool.close()
    pool.join()
    return results


def plotresults(dat, ind, msk, grd, pos, sig, tau, dec, geo, algo):
    ratio = getratio(grd, geo['mask']['step'])
    _dat = pack(dat)
    _ind = pack(ind)
    _pos = pack(pos)
    _sig = pack(sig)
    _tau = pack(tau)
    _dec = pack(dec)
    _dat = stretch(_dat.copy(), ratio)
    _msk = msk.copy()
    npix = _ind.shape[0]
    for m in range(npix):
        if _pos[m] > 0:
            if algo['method'] == 'lsqr':
                plotlsqr(_dat[m], _ind[m], _msk, _pos[m], _sig[m], _tau[m], _dec[m], ratio, m)
            logging.info("Saved: tmp/plot/plot-" + str(m) + ".png")


def plotlsqr(dat, ind, msk, pos, sig, tau, dec, ratio, id):
    import matplotlib.pyplot as plt 
    from scipy import ndimage   
    kernel = footprint(int(sig * ratio), tau, dec / np.sqrt(ratio))
    _msk = np.abs(msk - 1)
    _msk = convolve(_msk, kernel)
    tmp = ndimage.shift(_msk, -pos * ratio, order=1)[0:dat.size]
    plt.figure(figsize=(16, 3))
    plt.subplot(211)
    plt.title(str(id) + ' ' + 
        ' ind=' + str(ind) + 
        ' pos=' + str(pos * ratio) + 
        ' sig=' + str(sig * ratio) + 
        ' tau=' + str(tau) + 
        ' dec=' + str(dec / np.sqrt(ratio)))
    plt.grid('on')
    plt.plot(dat)
    plt.subplot(212)
    tmp = (tmp - tmp.min()) / (tmp.max() - tmp.min())
    plt.plot(tmp, 'r')
    stdmin = 2 * np.sqrt(dat.min())
    dat = dat  - (dat.min() + stdmin)
    stdmax = 2 * np.sqrt(dat.max())
    dat = dat / (dat.max() - stdmax)
    plt.plot(dat)
    plt.grid('on')
    plt.tight_layout()
    if not os.path.exists('tmp/plot'):
        os.makedirs('tmp/plot')
    plt.savefig('tmp/plot/plot-' + str(id) + '.png')
    plt.close()
        

def getratio(grid, step):
    ratio = step / (grid[1] - grid[0])
    logging.info("Ratio of signal-to-mask discretization is calculated: " + str(ratio))
    return ratio


def stretch(values, ratio):
    dx, dy = values.shape
    numimgs = int(dy * ratio)
    values = cv2.resize(values, (numimgs, dx))
    logging.info("Signal and mask discretizations are matched.")
    return values


def calibrate(data, ind, pos, sig, tau, dec, geo):
    _geo = [geo] * len(data)
    _dist = np.arange(*geo['mask']['calibrate']['dist'])
    dist = 1.165
    maximum = 0
    for k in range(_dist.size):
        _dis = [_dist[k]] * len(data)
        depth = runpar(_calibrate, [data, ind, pos, sig, tau, dec, _geo, _dis], len(data))
        saveplt('tmp/dep-dis/dep-' + str(k), depth, geo['source']['grid'])
        if np.max(sum(depth)) > maximum:
            maximum = np.max(sum(depth))
            dist = _dist[k]
    return dist


def _calibrate(args):
    # Unpack arguments
    dat = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    tau = args[4]
    dec = args[5]
    geo = args[6]
    geo['mask']['dist'] = args[7]

    # Source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')

    npixs = dat.shape[0]
    for m in range(npixs):
        # Calculate the signal footprint along the beam
        fp = _footprint(dat[m], ind[m], pos[m], sig[m], tau[m], dec[m], geo)
        depth += fp
    return depth


def footprintcor(geo, sig):
    p = geo['detector']['size'][0] / geo['detector']['shape'][0] * 1e3 # 200
    d = geo['mask']['dist'] # 1.17
    c = geo['detector']['pos'][2] # 513.140
    b = c - d # 511.969
    a = (sig * c - p * d) / (p - sig)
    return p * a / (a + c)


def resolve(dat, ind, pos, sig, tau, dec, geo):
    geo = [geo] * len(dat)
    results = runpar(_resolve, [dat, ind, pos, sig, tau, dec, geo], len(dat))
    depth = [None] * len(dat)
    laue = [None] * len(dat)
    for m in range(len(dat)):
        depth[m] = results[m][0]
        laue[m] = results[m][1]
    return depth, laue


def _resolve(args):
    # Unpack arguments
    dat = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    tau = args[4]
    dec = args[5]
    geo = args[6]
    
    # Number of pixels to be processed
    npixs = dat.shape[0]

    # Discrete grid along the source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')
    laue = np.zeros((npixs, len(grid)), dtype='float32')
    for m in range(npixs):
        # Calculate the signal footprint along the beam
        fp = _footprint(dat[m], ind[m], pos[m], sig[m], tau[m], dec[m], geo)
        laue[m] += fp
        depth += fp
    return depth, laue


def _footprint(dat, ind, pos, sig, tau, dec, geo):
    # Detector pixel position
    p0 = pix2pos(ind, geo)

    # Source beam
    s1 = np.array([-100, 0, 0], dtype='float32')
    s2 = np.array([100, 0, 0], dtype='float32')

    # Points on the mask for ray tracing
    p1 = np.array([-pos * 1e-3, geo['mask']['dist'], -100], dtype='float32')
    p2 = np.array([-pos * 1e-3, geo['mask']['dist'], 100], dtype='float32')
    intersection = intersect(s1, s2, p0, p1, p2)

    fp = __footprint(dat, sig, tau, dec, geo, intersection)
    return fp


def __footprint(dat, sig, tau, dec, geo, intersection):
    # Demagnification of the footprint
    _sig = footprintcor(geo, sig)

    # A fine discrete grid along the source beam
    gr = geo['source']['grid']
    grid = np.arange(gr[0], gr[1], gr[2])

    # Initialize fine footprint
    fp = np.zeros((len(grid), ), dtype='float32')

    # Index along the beam
    dx = np.argmin(np.abs(grid - intersection))
    if dx > 0 and dx < len(grid) - 1:
        fp[dx] = dat.max() - 2 * np.sqrt(dat.max())
        kernel = footprint(int(_sig * 1e-3 / gr[2]), tau, dec / (1e-3 / gr[2]))
        fp = convolve(fp, kernel)
    return fp


def intersect(s1, s2, p0, p1, p2):
    s12 = (s2 - s1) / np.linalg.norm(s2 - s1)
    p01 = (p1 - p0) / np.linalg.norm(p1 - p0) # vector A
    p02 = (p2 - p0) / np.linalg.norm(p2 - p0) # vector B
    vec = np.cross(p01, p02) # Normal vector
    t = np.dot(vec, s1 - p0) / np.dot(-s12, vec)
    intersect = s1[0] + t * s12[0]
    return intersect


def pix2pos(index, geo):
    """Returns the coordinate of a given detector pixel index."""
    dx = geo['detector']['shape'][0]
    dy = geo['detector']['shape'][1]    
    sx = geo['detector']['size'][0]
    sy = geo['detector']['size'][1]
    xp = (index[0] - 0.5 * (dx - 1)) * sx / dx
    yp = (index[1] - 0.5 * (dy - 1)) * sy / dy
    zp = 0
    xp += geo['detector']['pos'][0]
    yp += geo['detector']['pos'][1]
    zp += geo['detector']['pos'][2]
    xpoint = np.array([xp, yp, zp])
    rodrot = rotmatrix(geo['detector']['rot'])
    return np.dot(rodrot, xpoint)


def rotmatrix(rotation):
    """Returns the Rorigues's rotation matrix."""
    theta = np.sqrt(np.power(rotation, 2).sum())
    rx = rotation[0] / theta
    ry = rotation[1] / theta
    rz = rotation[2] / theta
    ctheta = np.cos(theta)
    stheta = np.sin(theta)
    c1 = 1 - ctheta
    matrix = np.zeros((3, 3), dtype='float32')
    matrix[0, 0] = ctheta + rx*rx*c1
    matrix[0, 1] = rx*ry*c1 - rz*stheta
    matrix[0, 2] = ry*stheta + rx*rz*c1
    matrix[1, 0] = rz*stheta + rx*ry*c1
    matrix[1, 1] = ctheta + ry*ry*c1
    matrix[1, 2] = -rx*stheta + ry*rz*c1
    matrix[2, 0] = -ry*stheta + rx*rz*c1
    matrix[2, 1] = rx*stheta + ry*rz*c1
    matrix[2, 2] = ctheta + rz*rz*c1
    return matrix


def maskcor(ind, mask, geo):
    correction = np.zeros(ind.shape, dtype='float32')
    npix = ind.shape[0]
    for m in range(npix):
        point = pix2pos(ind[m], geo)[0:2]
        orient = [0, geo['mask']['dist']]
        theta = np.arccos(
            np.dot(point, orient) / 
                (np.linalg.norm(point) * np.linalg.norm(orient)))
        displacement = np.tan(theta)
        correction[m] = (displacement * 
            mask['thickness'] / mask['resolution'])
    return np.floor(correction).astype('int32')
