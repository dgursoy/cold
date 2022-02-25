#!/usr/bin/env python3

import os
import numpy as np
import logging
from scipy import signal, ndimage, optimize, interpolate
from scipy.interpolate import BSpline
from cold import pack, partition, smooth, loadsingle, load
from skimage.feature import blob_log
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'





def decode(data, ind, mask, comp, geo, algo, debug=False):
    """Decodes the position and pixel footprints and their positions
    on the mask using coded measurement data."""

    # Initialize mask
    mask = invert(mask)
    logging.info('Initialized: Mask')

    # Initialize pos
    pos = np.zeros((data.shape[0], ), dtype='float32')
    logging.info('Initialized: Positions')

    # Initialize sig
    maxsize = algo['sig']['init']['maxsize']
    avgsize = algo['sig']['init']['avgsize']
    sig = np.zeros(maxsize, dtype='float32')
    first = int((maxsize - 1) * 0.5 - avgsize * 0.5)
    window = signal.windows.hann(int(avgsize))
    window /= sum(window)
    sig[first:first+avgsize] = window
    factor = 1 / geo['mask']['resolution']
    sig = stretcharr(sig, factor)
    sig = np.tile(sig, (data.shape[0], 1))
    logging.info('Initialized: Signal')

    # Initialize scl
    factor = geo['scanner']['step'] / geo['mask']['resolution']
    beta = (90 - geo['scanner']['angle']) * np.pi / 180 
    weighty = 1 / np.cos(geo['mask']['focus']['angley'] * np.pi / 180)
    weightz = np.cos(geo['mask']['focus']['anglez'] * np.pi / 180)
    sclscan = np.ones((ind.shape[0], ), dtype='float32')
    for m in range(ind.shape[0]):
        p0 = pix2pos(ind[m], geo) # [<->, dis2det, v^]
        alpha = np.arctan2(p0[0], p0[1]) #+25
        sclscan[m] = np.sin(beta) + np.cos(beta) * np.tan(alpha)
    scl = weighty * weightz * factor * geo['scanner']['step'] * sclscan
    logging.info('Initialized: Scales')

    # Initialize bases
    size = sig.shape[1]
    order = algo['sig']['order']
    scale = algo['sig']['scale']
    grid = np.arange(0, size + 1, dtype='int')
    knots = np.arange(0, size + 1, scale, dtype='int')
    ncoefs = knots.size - order - 1
    base = np.zeros((size, ncoefs))
    for m in range(ncoefs):
        t = knots[m:m + order + 2] # knots
        b = BSpline.basis_element(t)
        inc = knots[m]
        while knots[m + order + 1] > inc:
            base[inc, m] = b.integrate(grid[inc], grid[inc + 1])
            inc += 1
    logging.info('Initialized: Spline bases')

    # Partitioning
    ind = partition(ind, comp['workers'])
    pos = partition(pos, comp['workers'])
    sig = partition(sig, comp['workers'])
    scl = partition(scl, comp['workers'])
    data = partition(data, comp['workers'])
    datasize = data[0].shape[0] * data[0].shape[1] * 4e-6 # [MB]
    logging.info(
        "Data partitioned: " +
        "{} blocks of {}, {:.2f} MB each, {:.2f} GB total".format(
            comp['workers'], 
            data[0].shape, 
            datasize,
            datasize * comp['workers'] * 1e-3))
    logging.info('Partitioning completed')

    if comp['server'] == 'local':
        # Pack arguments as list and run   
        mask = [mask] * comp['workers']
        algo = [algo] * comp['workers']
        base = [base] * comp['workers']

        args = [data, mask, pos, sig, scl, algo, base]
        results = runpar(_decode, args, comp['workers'])

        # Unpack results and rescale them
        for m in range(comp['workers']):
            pos[m] = results[m][0]
            sig[m] = results[m][1]
        
    elif comp['server'] == 'remote':
        from funcx.sdk.client import FuncXClient
        import time

        fxc = FuncXClient()
        fxc.version_check()
        info = fxc.get_endpoint_status(comp['functionid'])
        if info['status'] != 'online': 
            raise RuntimeError("End point is not online!")
        fid = fxc.register_function(_decode, description=f"Decoder.")

        # Create a batch function tasks
        job = fxc.create_batch()
        for m in range(comp['workers']):
            args = [data[m], mask, pos[m], sig[m], scl[m], algo, base]
            job.add(args, endpoint_id=comp['functionid'], function_id=fid)

        # Run the batch function tasks
        batch_task_ids = fxc.batch_run(job)

        # Info
        counter = 0
        while True: 
            batch_task = fxc.get_batch_result(batch_task_ids)
            running_tasks = [s for s in batch_task 
                if batch_task[s]['status'] != 'success']
            finished_tasks = [s for s in batch_task 
                if batch_task[s]['status'] == 'success']
            if not running_tasks: break
            else: 
                logging.info(
                    f"Total remote runtime (sec): {counter}; " +
                    "Tasks are still running: " + 
                    f"{len(finished_tasks)} / {len(running_tasks)}")
            time.sleep(10) 
            counter += 10

        # Unpack results and rescale them
        for m, task in zip(range(comp['workers']), finished_tasks): 
            pos[m] = batch_task[task]['result'][0]
            sig[m] = batch_task[task]['result'][1]

    if debug is True:
        plotresults(data, ind, mask, pos, sig, scl, geo, algo)

    return pack(pos), pack(sig)


def _decode(args):
    data = args[0]
    mask = args[1]
    pos = args[2]
    sig = args[3]
    scl = args[4]
    algo = args[5]
    base = args[6]

    from cold import pixdecode
    import logging
    for m in range(data.shape[0]):
        pos[m], sig[m] = pixdecode(
            data[m], mask, pos[m], sig[m], scl[m], algo, base)
        logging.info('Pixel decoded: ' +
            str(m) + '/' + str(data.shape[0] - 1) + 
            ' pos=' + str(pos[m].squeeze()) + 
            ' scales=' + str(scl[m].squeeze()))
    return pos, sig


def pixdecode(data, mask, pos, sig, scale, algo, bases):
    """The main function for decoding pixel data."""
    data = normalize(data)
    data = ndimage.zoom(data, scale, order=1)
    pos = posrecon(data, mask, pos, sig, algo)
    sig = sigrecon(data, mask, pos, sig, algo, bases)
    return pos, sig


def posrecon(data, mask, pos, sig, algo):
    sim = signal.convolve(mask, sig, 'same')
    costsize = sim.size - data.size
    cost = np.zeros((costsize), dtype='float32')
    for m in range(costsize):
        if algo['pos']['method'] == 'lsqr':
            cost[m] = (np.sum(np.power(sim[m:m+data.size] - data, 2)) + 
                algo['pos']['regpar'] * np.sum(np.power(m - pos, 2)))
    try:
        pos = np.where(cost.min() == cost)[0][0]
    except IndexError:
        pass
    return pos


def sigrecon(data, mask, pos, sig, algo, base):
    first = int((sig.size - 1) / 2)
    last = int(mask.size + first - sig.size - data.size)
    if pos > first and pos < last:
        kernel = np.zeros((data.size, sig.size), dtype='float32')
        for m in range(data.size):
            begin = pos - first + m - 1
            end =  begin + sig.size
            kernel[m] = mask[begin:end]
        if algo['sig']['method'] == 'splines':
            coefs = optimize.nnls(np.dot(kernel, base), data)[0][::-1]
            sig = np.dot(base, coefs)
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


def smooth(img, sigma):
    img = ndimage.gaussian_filter(img, sigma)
    logging.info("Image smoothed.")
    return img


def blobsearch(img):
    blobs = blob_log(255. * img / img.max(), threshold=0.3)
    logging.info("Number of blobs found: " + str(blobs.shape[0]))
    return blobs


def sortblobs(blobs, img):
    nblobs = blobs.shape[0]
    arr = np.zeros((nblobs, ), dtype='float32')
    for m in range(nblobs):
        arr[m] = img[int(blobs[m, 0]), int(blobs[m, 1])]
    blobs = blobs[np.argsort(arr)]
    return blobs


def initializepos(file, mask, geo, algo, threshold=1000, debug=True):
    # Load an image
    img = loadsingle(file, id=file['range'][0])
    img[img < threshold] = 0
    img[img == 65535] = 0
    img = np.array(img, dtype='float32')
    img = signal.medfilt2d(img, kernel_size=3)
    
    # Find blob coordinates
    imgs = smooth(img.copy(), 2)
    blobs = blobsearch(imgs)
    blobs = sortblobs(blobs, imgs)

    # Find max points
    points = np.zeros(blobs[:, 0:2].shape, dtype='int')
    for m in range(blobs.shape[0]):
        x = int(blobs[m, 0])
        y = int(blobs[m, 1])
        a = img[x - 8: x + 8, y - 8: y + 8]
        i = np.unravel_index(np.argmax(a, axis=None), a.shape)
        if x >= 8:
            points[m, 0] = x - 8 + i[0] + 1
        if y >= 8:
            points[m, 1] = y - 8 + i[1] + 1

    # Find positions
    data, ind = load(file, collapsed=True, partitioned=True, index=points)
    algo['iter'] = 2
    algo['sig']['scale'] = 10
    algo['sig']['init']['avgsize'] = 40
    pos, sig, scales, tmp = decode(data, ind, mask, geo, algo)
    pos = pack(pos)
    ind = pack(ind)

    # Interpolate positions
    if algo['pos']['init'] == 'nearest':
        x = np.arange(0, 2048, 1)
        y = np.arange(0, 2048, 1)
        xx, yy = np.meshgrid(x, y)
        pos0 = interpolate.griddata(ind, pos, (yy, xx), method='nearest')
    elif algo['pos']['init'] == 'splines':
        xgrid = np.mgrid[0:2048:2048j, 0:2048:2048j]
        xflat = xgrid.reshape(2, -1).T
        # linear thin_plate_spline cubic quintic 
        RR = interpolate.RBFInterpolator(
            ind, pos, kernel='thin_plate_spline') 
        yflat = RR(xflat)
        pos0 = yflat.reshape(2048, 2048)

    # Testing
    if debug is True:
        import dxchange
        test = np.zeros(img.shape)
        for m in range(blobs.shape[0]):
            test[points[m, 0], points[m, 1]] = 1
        test = smooth(test, 2)
        dxchange.write_tiff(pos0 +  test / test.max() * 1e3, 'tmp/pos0.tiff')
        dxchange.write_tiff(pos0, 'tmp/pos0.tiff')
    return pos0


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


def stretcharr(arr, factor, order=1):
    arr = arr.copy()
    if len(arr.shape) == 2:
        zoom = [1, factor]
    elif len(arr.shape) == 1:
        zoom = factor
    arr = ndimage.zoom(arr, zoom, order=order)
    return arr / factor


def plotresults(dat, ind, msk, pos, sig, scl, geo, algo):
    _dat = pack(dat)
    _ind = pack(ind)
    _pos = pack(pos)
    _sig = pack(sig)
    _scl = pack(scl)
    _msk = msk.copy()
    npix = _ind.shape[0]
    for m in range(npix):
        if _pos[m] > 0:
            if algo['pos']['method'] == 'lsqr':
                plotlsqr(
                    _dat[m], _ind[m], _msk, 
                    _pos[m], _sig[m], _scl[m], m)
            logging.info("Saved: tmp/plot/plot-" + str(m) + ".png")


def plotlsqr(dat, ind, msk, pos, sig, scl, id):
    import matplotlib.pyplot as plt 
    plt.figure(figsize=(8, 4))
    
    dat = normalize(dat)
    dat = ndimage.zoom(dat, scl, order=1)
    _sig = sig / sig.max()

    msk = invert(msk)
    _dat = np.zeros(msk.shape)
    _dat[int(pos):int(pos+dat.size)] = dat

    plt.subplot(311)
    plt.title(str(id) + ' ' + ' ind=' + str(ind) + ' pos=' + str(pos))
    plt.step(_sig, 'darkorange')
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    msk = signal.convolve(msk, sig, 'same')

    plt.subplot(312)
    plt.step(msk[300:6060], 'tab:blue')
    plt.step(_dat[300:6060], 'tab:red')
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    _msk = ndimage.shift(msk, -pos, order=1)[0:dat.size]

    plt.subplot(313)
    plt.step(dat, 'tab:red')
    plt.step(_msk, 'tab:blue')
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    plt.tight_layout()
    if not os.path.exists('tmp/plot'):
        os.makedirs('tmp/plot')
    plt.savefig('tmp/plot/plot-' + str(id) + '.png', dpi=240)
    plt.close()


# def calibrate(data, ind, pos, sig, shift, geo):
#     # Number pf processes
#     chunks = len(data)

#     # Initialize 
#     dist = np.arange(*geo['mask']['calibrate']['dist'])

#     maximum = 0
#     optdist = 0
#     for k in range(len(dist)):
#         # Pack arguments as list and run
#         geo['mask']['focus']['dist'] = dist[k]
#         args = packcalibrateargs(data, ind, pos, sig, shift, geo, chunks)
#         depth = runpar(_calibrate, args, chunks)

#         # Save results
#         saveplt('tmp/dep-dist/dep-' + str(k), depth, geo['source']['grid'])
#         if np.sqrt(np.sum(np.power(depth, 2))) > maximum:
#             maximum = np.sqrt(np.sum(np.power(depth, 2)))
#             optdist = dist[k]
#     logging.info("Optimum distance: " + str(optdist))
#     return optdist


# def packcalibrateargs(data, ind, pos, sig, shift, geo, chunks):
#     geo = [geo] * chunks
#     shift = [shift] * chunks
#     return [data, ind, pos, sig, shift, geo]


# def unpackcalibrateargs(args):
#     data = args[0]
#     ind = args[1]
#     pos = args[2]
#     sig = args[3]
#     shift = args[4]
#     geo = args[5]
#     return data, ind, pos, sig, shift, geo


# def _calibrate(args):
#     # Unpack arguments
#     data, ind, pos, sig, shift, geo = unpackcalibrateargs(args)

#     # Source beam
#     grid = np.arange(*geo['source']['grid'])

#     # Initializations
#     depth = np.zeros((len(grid), ), dtype='float32')

#     # Number of pixels to be processed
#     npix = data.shape[0]

#     for m in range(npix):
#         # Calculate the signal footprint along the beam
#         fp = _footprint(data[m], ind[m], pos[m], sig[m], shift, geo)
#         depth += fp
#     return depth


def resolve(data, ind, pos, sig, geo, shift=0):
    """Resolves depth information."""
    # Number pf processes
    chunks = len(data)

    # Pack arguments as list and run
    shift = [shift] * chunks
    geo = [geo] * chunks
    args = [data, ind, pos, sig, shift, geo]
    results = runpar(_resolve, args, chunks)

    # Unpack results
    depth = [None] * chunks
    laue = [None] * chunks
    for m in range(chunks):
        depth[m] = results[m][0]
        laue[m] = results[m][1]
    return depth, laue


def _resolve(args):
    # Unpack arguments
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    shift = args[4]
    geo = args[5]
    
    # Number of pixels to be processed
    npix = data.shape[0]

    # Discrete grid along the source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')
    laue = np.zeros((npix, len(grid)), dtype='float32')
    for m in range(npix):
        # Calculate the signal footprint along the beam
        fp = _footprint(data[m], ind[m], pos[m], sig[m], shift, geo)
        laue[m] += fp
        depth += fp
    return depth, laue


def _rotpoint2d(x, y, cx, cy, angle):
    px = x - cx
    py = y - cy
    xnew = px * np.cos(angle) - py * np.sin(angle)
    ynew = px * np.sin(angle) + py * np.cos(angle)
    px = xnew + cx
    py = ynew + cy
    return px, py


def rotpoint2d(x, y, angle):
    xnew = x * np.cos(angle) - y * np.sin(angle)
    ynew = x * np.sin(angle) + y * np.cos(angle)
    return xnew, ynew


def _footprint(dat, ind, pos, sig, shift, geo):
    # Detector pixel position
    p0 = pix2pos(ind, geo) # [<->, dis2det, v^]

    # Source beam
    s1 = np.array([-100, 
        -geo['source']['offset'] * 0.5, 0], dtype='float32')
    s2 = np.array([100, 
        geo['source']['offset'] * 0.5, 0], dtype='float32')

    # Scaling of positions and signals
    pos *= geo['mask']['resolution']
    sig = ndimage.zoom(sig, geo['mask']['resolution'], order=1)

    # Mask position
    xx = geo['mask']['focus']['cenx'] + shift
    yy = geo['mask']['focus']['dist']
    zz = geo['mask']['focus']['cenz']

    # Points on the mask for ray tracing
    p1 = np.array([-pos * 1e-3 + shift, yy, 100], dtype='float32')
    p2 = np.array([-pos * 1e-3 + shift, yy, -100], dtype='float32')

    # Rotate a rigid body
    anglex = geo['mask']['focus']['anglex'] * np.pi / 180
    angley = geo['mask']['focus']['angley'] * np.pi / 180
    anglez = geo['mask']['focus']['anglez'] * np.pi / 180
    p1[1], p1[2] = _rotpoint2d(p1[1], p1[2], yy, zz, anglex)
    p2[1], p2[2] = _rotpoint2d(p2[1], p2[2], yy, zz, anglex)
    p1[0], p1[2] = _rotpoint2d(p1[0], p1[2], xx, zz, angley)
    p2[0], p2[2] = _rotpoint2d(p2[0], p2[2], xx, zz, angley)
    p1[0], p1[1] = _rotpoint2d(p1[0], p1[1], xx, yy, anglez)
    p2[0], p2[1] = _rotpoint2d(p2[0], p2[1], xx, yy, anglez)

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
    d = geo['mask']['focus']['dist'] # 1.17
    c = geo['detector']['pos'][2] # 513.140
    b = c - d # 511.969
    try:
        a = (sig.size * c - p * d) / (p - sig.size)
        factor =  p * a / (a + c) / sig.size  / (1e3 * gr[2])
    except ZeroDivisionError:
        factor = 1
    return stretcharr(sig, factor, order=1)


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
        orient = [0, geo['mask']['focus']['dist']]
        theta = np.arccos(
            np.dot(point, orient) / 
                (np.linalg.norm(point) * np.linalg.norm(orient)))
        displacement = np.tan(theta)
        correction[m] = (displacement * 
            mask['thickness'] / mask['resolution'])
    return np.floor(correction).astype('int32')
