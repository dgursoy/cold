#!/usr/bin/env python3

import numpy as np
import logging
from scipy import signal, ndimage, optimize, linalg
from cold import pack, saveplt, partition
import multiprocessing
import os


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'



def decode(data, mask, geo, algo):
    """Decodes the position and pixel footprints and their positions
    on the mask using coded measurement data."""

    # Number pf processes
    chunks = len(data)

    # Number of total pixels
    npix = 0
    for item in data:
        npix += item.shape[0]

    # Scaling factor
    factor = stretchfactor(
        geo['mask']['resolution'], 
        geo['mask']['step'])

    # Initialize parameters
    data = initdata(data, factor)
    mask = initmask(mask)
    pos = initpos(algo, chunks, npix, factor)
    sig = initsig(algo, chunks, npix, factor)

    # Pack arguments as list and run
    args = packdecodeargs(data, mask, pos, sig, algo, chunks)
    results = runpar(_decode, args, chunks)
    
    # Unpack results and rescale them
    for m in range(chunks):
        pos[m] = results[m][0] / factor
        sig[m] = stretcharr(results[m][1], 1 / factor)
    return pos, sig


def stretchlist(data, factor):
    data = data.copy()
    chunks = len(data)
    for m in range(chunks):
        data[m] = stretcharr(data[m], factor)
    return data


def stretcharr(arr, factor, order=1):
    arr = arr.copy()
    if len(arr.shape) == 2:
        zoom = [1, factor]
    elif len(arr.shape) == 1:
        zoom = factor
    arr = ndimage.zoom(arr, zoom, order=order)
    return arr
        
    
def initdata(data, factor):
    return stretchlist(data, factor)


def initmask(mask):
    return invert(mask)


def initpos(algo, chunks, npix, factor):
    pos = np.array(algo['pos']['init'] * factor, dtype='float32')
    pos = np.zeros((npix, ), dtype='float32')
    return partition(pos, chunks)


def initsig(algo, chunks, npix, factor):
    maxsize = algo['sig']['init']['maxsize']
    avgsize = algo['sig']['init']['avgsize']
    sig = np.zeros(maxsize, dtype='float32')
    first = int((maxsize - 1) * 0.5 - avgsize * 0.5)
    sig[first:first+avgsize] = footprnt(avgsize)
    sig = stretcharr(sig, factor) / factor
    sig = np.tile(sig, (npix, 1))
    return partition(sig, chunks)


def packdecodeargs(data, mask, pos, sig, algo, chunks):
    mask = [mask] * chunks
    algo = [algo] * chunks
    logging.info("Packed arguments as list")
    return [data, mask, pos, sig, algo]


def unpackdecodeargs(args):
    data = args[0]
    mask = args[1]
    pos = args[2]
    sig = args[3]
    algo = args[4]
    logging.info("Unpacked list.")
    return data, mask, pos, sig, algo


def stretchfactor(resolution, step):
    factor = step / resolution
    logging.info("Factor is calculated: " + str(factor))
    return factor


def _decode(args):
    data, mask, pos, sig, algo = unpackdecodeargs(args)
    npix = data.shape[0]
    for m in range(npix):
        pos[m], sig[m] = pixdecode(data[m], mask, pos[m], sig[m], algo)
        logging.info('Pixel decoded: ' +
            str(m) + '/' + str(npix - 1) + 
            ' pos=' + str(pos[m].squeeze()))
    return pos, sig


def pixdecode(data, mask, pos, sig, algo):
    """The main function for decoding pixel data."""
    data = normalize(data)
    pos = posrecon(data, mask, pos, sig, algo)
    sig = sigrecon(data, mask, pos, sig, algo)
    return pos, sig


def posrecon(data, mask, pos, sig, algo):
    sim = signal.convolve(mask, sig, 'same')
    costsize = sim.size - data.size
    cost = np.zeros((costsize), dtype='float32')
    for m in range(costsize):
        if algo['pos']['method'] == 'lsqr':
            cost[m] = np.sum(np.power(sim[m:m+data.size] - data, 2))
        pos = np.where(cost.min() == cost)[0][0]
    return pos


def sigrecon(data, mask, pos, sig, algo):
    first = int((sig.size - 1) / 2)
    last = int(mask.size + first - sig.size - data.size)
    if pos > first and pos < last:
        kernel = np.zeros((data.size, sig.size), dtype='float32')
        for m in range(data.size):
            begin = pos - first + m
            end =  begin + sig.size
            kernel[m] = mask[begin:end]
        if algo['sig']['method'] == 'nnls':
            sig = optimize.nnls(kernel, data)[0][::-1]
        if algo['sig']['method'] == 'pinv':
            ikernel = linalg.pinv(kernel, algo['sig']['init']['atol'])
            sig = np.dot(ikernel, data)[::-1]
        sig /= sig.sum()
    else:
        sig *= 0
    return sig


def normalize(data):
    data = data.copy()
    diff = max(data) - min(data)
    stdmin = 2 * np.sqrt(min(data))
    if stdmin > diff * 0.5:
        stdmin = 0
    data = data - (min(data) + stdmin)
    stdmax = 2 * np.sqrt(max(data))
    data = data / (max(data) - stdmax)
    return data


def footprnt(val):
    """Create a discrete signal from its parameters."""
    window = signal.windows.hann(int(val))
    window /= sum(window)
    return window


def invert(mask):
    return np.abs(mask.copy() - 1)


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


def plotresults(dat, ind, msk, pos, sig, geo, algo):
    factor = stretchfactor(
        geo['mask']['resolution'], 
        geo['mask']['step'])
    _dat = pack(dat)
    _ind = pack(ind)
    _pos = pack(pos)
    _sig = pack(sig)
    _dat = stretcharr(_dat.copy(), factor)
    _msk = msk.copy()
    npix = _ind.shape[0]
    for m in range(npix):
        if _pos[m] > 0:
            if algo['pos']['method'] == 'lsqr':
                plotlsqr(_dat[m], _ind[m], _msk, _pos[m], _sig[m], factor, m)
            logging.info("Saved: tmp/plot/plot-" + str(m) + ".png")


def plotlsqr(dat, ind, msk, pos, sig, factor, id):
    import matplotlib.pyplot as plt 
    from scipy import ndimage   
    sig = ndimage.zoom(sig, factor)
    _msk = invert(msk)
    _msk = signal.convolve(_msk, sig, 'same')
    tmp = ndimage.shift(_msk, -pos * factor, order=1)[0:dat.size]
    plt.figure(figsize=(16, 3))
    plt.subplot(211)
    plt.title(str(id) + ' ' + ' ind=' + str(ind))
    plt.grid('on')
    plt.plot(dat)
    plt.subplot(212)
    tmp = (tmp - tmp.min()) / (tmp.max() - tmp.min())
    plt.plot(tmp, 'r')
    dat = normalize(dat)
    plt.plot(dat)
    plt.grid('on')
    plt.tight_layout()
    if not os.path.exists('tmp/plot'):
        os.makedirs('tmp/plot')
    plt.savefig('tmp/plot/plot-' + str(id) + '.png')
    plt.close()


def calibratetilt(data, ind, pos, sig, geo):
    # Number pf processes
    chunks = len(data)

    # Initialize 
    tiltx = np.arange(*geo['mask']['calibrate']['tiltx'])
    tilty = np.arange(*geo['mask']['calibrate']['tilty'])

    cost = np.zeros((len(tilty), len(tiltx)))
    maximum = 0
    opttiltx = 0
    opttilty = 0
    for m in range(len(tilty)):
        # Pack arguments as list and run
        geo['mask']['tilty'] = tilty[m]
        for n in range(len(tiltx)):
            geo['mask']['tiltx'] = tiltx[n]
            args = packcalibrateargs(data, ind, pos, sig, geo, chunks)
            depth = runpar(_calibrate, args, chunks)
            cost[m, n] = np.max(depth)
            # cost[k, m] = np.sqrt(np.sum(np.power(depth, 2)))

            # Save results
            saveplt('tmp/dep-tilt/dep-' + str(m) + '-' + str(n), depth, geo['source']['grid'])
            if cost[m, n] > maximum:
                maximum = cost[m, n]
                opttiltx = tiltx[n]
                opttilty = tilty[m]
            print (cost[m, n], maximum, opttiltx, opttilty)
    return opttiltx, opttilty


def calibratedist(data, ind, pos, sig, geo):
    # Number pf processes
    chunks = len(data)

    # Initialize 
    dist = np.arange(*geo['mask']['calibrate']['dist'])

    maximum = 0
    optdist = 0
    for k in range(len(dist)):
        # Pack arguments as list and run
        geo['mask']['dist'] = dist[k]
        args = packcalibrateargs(data, ind, pos, sig, geo, chunks)
        depth = runpar(_calibrate, args, chunks)

        # Save results
        saveplt('tmp/dep-dist/dep-' + str(k), depth, geo['source']['grid'])
        if np.sqrt(np.sum(np.power(depth, 2))) > maximum:
            maximum = np.sqrt(np.sum(np.power(depth, 2)))
            optdist = dist[k]
    return optdist


def packcalibrateargs(data, ind, pos, sig, geo, chunks):
    geo = [geo] * chunks
    return [data, ind, pos, sig, geo]


def unpackcalibrateargs(args):
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    geo = args[4]
    return data, ind, pos, sig, geo


def _calibrate(args):
    # Unpack arguments
    data, ind, pos, sig, geo = unpackcalibrateargs(args)

    # Source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')

    # Number of pixels to be processed
    npix = data.shape[0]

    for m in range(npix):
        # Calculate the signal footprint along the beam
        fp = _footprint(data[m], ind[m], pos[m], sig[m], geo)
        depth += fp
    return depth


def resolve(data, ind, pos, sig, geo):
    """Resolves depth information."""
    # Number pf processes
    chunks = len(data)

    # Pack arguments as list and run
    args = packresolveargs(data, ind, pos, sig, geo, chunks)
    results = runpar(_resolve, args, chunks)

    # Unpack results
    depth = [None] * chunks
    laue = [None] * chunks
    for m in range(chunks):
        depth[m] = results[m][0]
        laue[m] = results[m][1]
    return depth, laue


def packresolveargs(data, ind, pos, sig, geo, chunks):
    geo = [geo] * chunks
    return [data, ind, pos, sig, geo]


def unpackresolveargs(args):
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    geo = args[4]
    return data, ind, pos, sig, geo


def _resolve(args):
    # Unpack arguments
    data, ind, pos, sig, geo = unpackresolveargs(args)
    
    # Number of pixels to be processed
    npix = data.shape[0]

    # Discrete grid along the source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')
    laue = np.zeros((npix, len(grid)), dtype='float32')
    for m in range(npix):
        # Calculate the signal footprint along the beam
        fp = _footprint(data[m], ind[m], pos[m], sig[m], geo)
        laue[m] += fp
        depth += fp
    return depth, laue


def _footprint(dat, ind, pos, sig, geo):
    # Detector pixel position
    p0 = pix2pos(ind, geo) # [<->, dis2det, v^]

    # Source beam
    s1 = np.array([-100, 0, 0], dtype='float32')
    s2 = np.array([100, 0, 0], dtype='float32')

    # Points on the mask for ray tracing
    tx = geo['mask']['tiltx'] * (p0[0] - geo['mask']['cenx']) / geo['detector']['size'][0] * 0.5
    ty = geo['mask']['tilty'] * (p0[2] - geo['mask']['ceny']) / geo['detector']['size'][1] * 0.5
    # print (pos, tx, p0[0], p0[1], p0[2])
    p1 = np.array([-(pos + tx + ty) * 1e-3, geo['mask']['dist'], -100], dtype='float32')
    p2 = np.array([-(pos + tx + ty) * 1e-3, geo['mask']['dist'], 100], dtype='float32')
    intersection = intersect(s1, s2, p0, p1, p2)

    # Discrete grid along the source beam
    gr = geo['source']['grid']
    grid = np.arange(gr[0], gr[1], gr[2])
    extgrid = np.arange(gr[0] - 1.0, gr[1] + 1.0, gr[2])

    # Demagnification of the footprint
    _sig = footprintcor(geo, sig)

    # Initialize fine footprint
    fp = np.zeros((len(extgrid), ), dtype='float32')

    # Index along the beam
    dx = np.argmin(np.abs(extgrid - intersection))
    first = int((_sig.size - 1) / 2)
    if extgrid[dx] > gr[0] and extgrid[dx] < gr[1]:
        _fp = _sig * (dat.max() - 2 * np.sqrt(dat.max()))
        begin = dx - first
        end = begin + _sig.size
        fp[begin:end] = _fp

    x = np.argmin(np.abs(extgrid - gr[0]))
    return fp[x:x+grid.size]


def footprintcor(geo, sig):
    gr = geo['source']['grid']
    p = geo['detector']['size'][0] / geo['detector']['shape'][0] * 1e3 # 200
    d = geo['mask']['dist'] # 1.17
    c = geo['detector']['pos'][2] # 513.140
    b = c - d # 511.969
    a = (sig.size * c - p * d) / (p - sig.size)
    factor =  p * a / (a + c) / sig.size  / (1e3 * gr[2])
    return stretcharr(sig, factor, order=3) / factor


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
