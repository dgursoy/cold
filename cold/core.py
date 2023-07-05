#!/usr/bin/env python3

import os
import numpy as np
import logging
from scipy import signal, ndimage, optimize, interpolate, linalg
from scipy.interpolate import BSpline
from cold import pack, partition, loadsingle, load, saveplt
from cold import mask as cmask
from skimage.feature import blob_log
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'





def decode(data, ind, comp, geo, algo, pos=None, debug=False):
    """Decodes the position and pixel footprints and their positions
    on the mask using coded measurement data."""

    # Initialize pos
    if pos is None:
        pos = np.zeros((data.shape[0], ), dtype='float32')
        logging.info('Initialized: Positions')

    # Initialize ene
    ene = np.zeros((data.shape[0], ), dtype='float32')
    logging.info('Initialized: Energy')

    pathlen = np.zeros((data.shape[0], ), dtype='float32')
    logging.info('Initialized: Pathlen')
    
    # Initialize new scl
    scl = np.ones((ind.shape[0], ), dtype='float32')
    factor = (geo['scanner']['step'] /
        geo['mask']['resolution'] * geo['mask']['stretch'])

    # # Rotation vector (intrinsic-zyx)
    # alpha = geo['mask']['focus']['anglez'] * np.pi / 180
    # beta = geo['mask']['focus']['angley'] * np.pi / 180
    # gamma = geo['mask']['focus']['anglex'] * np.pi / 180
    # rotmat = np.zeros((3, 3), dtype='float32')
    # rotmat[0, 0] = np.cos(alpha) * np.cos(beta)
    # rotmat[0, 1] = np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)
    # rotmat[0, 2] = np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma)
    # rotmat[1, 0] = np.sin(alpha) * np.cos(beta)
    # rotmat[1, 1] = np.sin(alpha) * np.sin(beta) * np.sin(gamma) - np.cos(alpha) * np.cos(gamma)
    # rotmat[1, 2] = np.sin(alpha) * np.sin(beta) * np.cos(gamma) + np.cos(alpha) * np.sin(gamma)
    # rotmat[2, 0] = -np.sin(beta)
    # rotmat[2, 1] = np.cos(beta) * np.sin(gamma)
    # rotmat[2, 2] = np.cos(beta) * np.cos(gamma)

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

    # Mask center step 1
    mc1 = np.array([
        geo['mask']['focus']['cenx'], 
        geo['mask']['focus']['dist'],
        geo['mask']['focus']['cenz']], dtype='float32')

    # Scan direction
    rodrot = rotmatrix(geo['scanner']['rot'])
    direction = np.dot(rodrot, geo['scanner']['axis'])
    mc1 = mc1 + direction * geo['mask']['shift']
    mc2 = mc1 + direction * geo['scanner']['step']

    # # Sample origin
    s0 = np.array([geo['source']['offset'], 0, 0], dtype='float32')

    # Rotation of mask axes
    mx = np.array([1, 0, 0], dtype='float32')
    my = np.array([0, 1, 0], dtype='float32')
    mz = np.array([0, 0, 1], dtype='float32')
    mx = np.dot(rotmat, mx)
    my = np.dot(rotmat, my)
    mz = np.dot(rotmat, mz) 

    for m in range(ind.shape[0]):
        # Detector pixel position
        p0 = pix2pos(ind[m], geo) # [<->, dis2det, v^]

        # Intersection of step 1
        i1 = intersect2(p0, s0, mc1, mc1 + mx, mc1 + mz)

        # Intersection of step 2
        i2 = intersect2(p0, s0, mc2, mc2 + mx, mc2 + mz)

        # Vector from i1 to i2
        vec = i1 - i2 + direction

        # Projection onto mask axis
        step = np.dot(vec, mx) 
        scl[m] = factor * step
    logging.info('Initialized: Scales')

    # Initialize sig
    maxsize = algo['sig']['init']['maxsize']
    sig = np.zeros(maxsize, dtype='float32')
    sig = np.tile(sig, (data.shape[0], 1))
    for m in range(ind.shape[0]):
        # Detector pixel position
        p0 = pix2pos(ind[m], geo) # [<->, dis2det, v^]

        angpix = np.arctan(p0[0] / p0[1]) 
        angmsk = geo['mask']['focus']['anglez'] * np.pi / 180.
        xxx = np.cos(angpix) / np.cos(angpix + angmsk)

        zoomin = 10
        _sig = np.zeros(maxsize * zoomin, dtype='float32')
        avgsize = algo['sig']['init']['avgsize'] * xxx * zoomin
        first = int((maxsize * zoomin - 1) * 0.5 - avgsize * 0.5)
        window = signal.windows.tukey(int(avgsize), alpha=0.1)
        window /= sum(window)
        _sig[first:first+int(avgsize)] = window
        sig[m] = ndimage.zoom(_sig * zoomin, 1 / zoomin, order=1)
        sig[m] /= sum(sig[m])

    factor = int(1 / geo['mask']['resolution'])
    sig = stretcharr(sig, factor)
    logging.info('Initialized: Signal')
    
    # # Initialize sig
    # maxsize = algo['sig']['init']['maxsize']
    # avgsize = algo['sig']['init']['avgsize']
    # sig = np.zeros(maxsize, dtype='float32')
    # first = int((maxsize - 1) * 0.5 - avgsize * 0.5)
    # window = signal.windows.tukey(int(avgsize), alpha=0.0)
    # window /= sum(window)
    # sig[first:first+avgsize] = window
    # factor = 1 / geo['mask']['resolution']
    # sig = stretcharr(sig, factor)
    # sig = np.tile(sig, (data.shape[0], 1))
    # logging.info('Initialized: Signal')
        
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
    if comp['server'] != 'proc':
        pos = partition(pos, comp['workers'])
        sig = partition(sig, comp['workers'])
        scl = partition(scl, comp['workers'])
        data = partition(data, comp['workers'])
        ind = partition(ind, comp['workers'])
        ene = partition(ene, comp['workers'])
        pathlen = partition(pathlen, comp['workers'])
        datasize = data[0].shape[0] * data[0].shape[1] * 4e-6 # [MB]
        logging.info(
            "Data partitioned: " +
            "{} blocks of {}, {:.2f} MB each, {:.2f} GB total".format(
                comp['workers'], 
                data[0].shape, 
                datasize,
                datasize * comp['workers'] * 1e-3))
        logging.info('Partitioning completed')

        # Pack arguments as list and run   
        algo = [algo] * comp['workers']
        base = [base] * comp['workers']
        geo = [geo] * comp['workers']

    args = [data, pos, sig, scl, algo, base, geo, ind, ene, pathlen]


    if comp['server'] == 'local':

        results = runpar(_decode, args, comp['workers'])

        # Unpack results and rescale them
        for m in range(comp['workers']):
            pos[m] = results[m][0]
            sig[m] = results[m][1]
            ene[m] = results[m][2]
            pathlen[m] = results[m][3]
        
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
            args = [data[m], pos[m], sig[m], scl[m], algo, base, geo, ind]
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
    
    elif comp['server'] == 'proc':
        results = _decode(args)
    
    if comp['server'] != 'proc':
        data = pack(data)
        pos = pack(pos)
        sig = pack(sig) 
        scl = pack(scl) 
        ind = pack(ind) 
        ene = pack(ene)
        pathlen = pack(pathlen)

    if debug == True:
        plotresults(data, ind, geo[0], pos, sig, scl, algo[0], ene)

    return pos, sig, scl, ene, pathlen


def _decode(args):
    data = args[0]
    pos = args[1]
    sig = args[2]
    scl = args[3]
    algo = args[4]
    base = args[5]
    geo = args[6]
    ind = args[7]
    ene = args[8]
    pathlen = args[9]

    from cold import pixdecode
    import logging
    for m in range(data.shape[0]):
        pos[m], sig[m], ene[m], pathlen[m] = pixdecode(
            data[m], pos[m], sig[m], scl[m], algo, base, geo, ind[m], ene[m], m)
        logging.info('Pixel decoded: ' +
            str(m) + '/' + str(data.shape[0] - 1) + 
            ' pos=' + str(pos[m].squeeze()) + 
            ' ene=' + str(ene[m].squeeze()) + 
            ' scales=' + str(scl[m].squeeze()))
    return pos, sig, ene, pathlen


def pixdecode(data, pos, sig, scale, algo, bases, geo, ind, ene, ix):
    """The main function for decoding pixel data."""
    msk = cmask.discmask(geo, ind)
    data = normalize(data)
    data = ndimage.zoom(data, scale, order=1)
    pos = posrecon(data, msk, pos, sig, algo)

    if algo['ene']['recon'] is True:
        ene, pathlen = enerecon(data, pos, sig, scale, algo, geo, ind)
        msk = cmask.discmask(geo, ind, exact=True, energy=ene)

    if algo['sig']['recon'] is True:
        sig = sigrecon(data, msk, pos, sig, algo, bases, ix)

    return pos, sig, ene, pathlen


def enerecon(data, pos, sig, scl, algo, geo, ind):
    _data = ndimage.zoom(data, scl, order=1)
    ex = algo['ene']['range']
    enerange = np.arange(*ex)
    costsize = len(enerange)
    cost = np.zeros((costsize), dtype='float32')
    pathlens = []
    for m in range(costsize):
        msk, pathlen = cmask.discmask(geo, ind, exact=True, energy=enerange[m], return_pathlen=True)
        pathlens.append(pathlen)
        sim = signal.convolve(msk, sig, 'same')
        _sim = ndimage.shift(sim, -pos, order=1)[0:_data.size]
        if algo['ene']['method'] == 'lsqr':
            cost[m] = np.sum(np.power(_data - _sim, 2))
    try:
        ii = np.where(cost.min() == cost)[0][0]
    except IndexError:
        ii = 0
        pass
    return enerange[ii], pathlens[ii]


def posrecon(data, msk, pos, sig, algo):
    sim = signal.convolve(msk, sig, 'same')
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


def sigrecon(data, msk, pos, sig, algo, base, ix):
    first = int((sig.size - 1) / 2)
    last = int(msk.size + first - sig.size - data.size)
    if pos > first and pos < last:
        kernel = np.zeros((data.size, sig.size), dtype='float32')
        for m in range(data.size):
            begin = pos - first + m - 1
            end =  begin + sig.size
            kernel[m] = msk[begin:end]
        if algo['sig']['method'] == 'splines':
            # np.save('tmp/analysis/kernel.npy', kernel)
            coefs = optimize.nnls(np.dot(kernel, base), data)[0][::-1]
            # np.save('tmp/analysis/spline.npy', np.dot(kernel, base))
            # np.save('tmp/analysis/base.npy', base)
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


def plotresults(dat, ind, geo, pos, sig, scl, algo, ene):
    npix = ind.shape[0]
    for m in range(npix):
        if algo['ene']['recon'] is True:
            msk = cmask.discmask(geo, ind[m], exact=True, energy=ene[m])
        else:
            msk = cmask.discmask(geo, ind[m])
        if pos[m] > 0:
            if algo['pos']['method'] == 'lsqr':
                plotlsqr(
                    dat[m], ind[m], msk, 
                    pos[m], sig[m], scl[m], m)
            logging.info("Saved: tmp/plot/plot-" + str(m) + ".png")


def plotlsqr(dat, ind, msk, pos, sig, scl, id):
    import matplotlib.pyplot as plt 
    import matplotlib 
    matplotlib.rcParams.update({
        "text.usetex": True, 
        'font.size': 11, 
        'font.family' : 'serif',
        })

    plt.figure(figsize=(6, 3))
    
    dat = normalize(dat)
    dat = ndimage.zoom(dat, scl, order=1)
    _sig = sig / sig.max()

    _dat = np.zeros(msk.shape)
    _dat[int(pos):int(pos+dat.size)] = dat

    plt.subplot(311)
    plt.title(str(id) + ' ' + ' ind=' + str(ind) + ' pos=' + str(pos))
    plt.step(_sig, 'orangered', linewidth=1.3)
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    sim = signal.convolve(msk.copy(), sig, 'same')
    sim /= sim.max()

    plt.subplot(312)
    plt.step(msk[300:6060], 'tab:green', linewidth=0.6)
    plt.step(sim[300:6060], 'tab:blue', linewidth=1)
    plt.step(_dat[300:6060], 'tab:red', linewidth=0.6)
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    _msk = ndimage.shift(msk, -pos, order=1)[0:dat.size]
    _sim = ndimage.shift(sim, -pos, order=1)[0:dat.size]

    plt.subplot(313)
    plt.step(_msk, 'tab:green', linewidth=0.6)
    plt.step(_sim, 'tab:blue', linewidth=1)
    plt.step(dat, 'tab:red', linewidth=0.6)
    plt.grid('on')
    plt.ylim((-0.5, 1.5))

    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.9, hspace=0.4)
    if not os.path.exists('tmp/plot'):
        os.makedirs('tmp/plot')
    filename = 'tmp/plot/plot-' + str(id) + '.png'
    incr = 0
    while os.path.exists(filename):
        filename = 'tmp/plot/plot-' + str(id + incr) + '.png'
        incr += 1
    plt.savefig(filename, dpi=480)
    plt.close()

def calibrate(data, ind, pos, sig, geo, comp, shift=0, depw=None):
    # Initialize 
    dist = np.arange(*geo['mask']['calibrate']['dist'])

    maximum = 0
    optdist = 0
    for k in range(len(dist)):
        geo['mask']['focus']['dist'] = dist[k]
        depth, laue = resolve(data, ind, pos, sig, geo, comp, shift=shift)

        # Save results
        saveplt('tmp/dep-dist/dep-' + str(k), depth, geo['source']['grid'], depw)
        if np.sqrt(np.sum(np.power(depth, 2))) > maximum:
            maximum = np.sqrt(np.sum(np.power(depth, 2)))
            optdist = dist[k]
    logging.info("Optimum distance: " + str(optdist))
    return optdist


def resolve(data, ind, pos, sig, geo, comp):
    """Resolves depth information."""

    # Pack arguments as list and run
    pos = partition(pos, comp['workers'])
    sig = partition(sig, comp['workers'])
    data = partition(data, comp['workers'])
    ind = partition(ind, comp['workers'])
    geo = [geo] * comp['workers']
    args = [data, ind, pos, sig, geo]
    results = runpar(_resolve, args, comp['workers'])

    # Unpack results
    depth = [None] * comp['workers']
    laue = [None] * comp['workers']
    for m in range(comp['workers']):
        depth[m] = results[m][0]
        laue[m] = results[m][1]
    depth = sum(depth)
    laue = pack(laue)
    return depth, laue


def _resolve(args):
    # Unpack arguments
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    geo = args[4]
    
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
    s1 = np.array([-1, -0.0068401066731 * 0.5, 0], dtype='float32')
    s2 = np.array([1, 0.0068401066731 * 0.5, 0], dtype='float32')

    # Scaling of positions and signals
    pos *= geo['mask']['resolution'] * 1e-3
    mask = cmask.discmask(geo, ind)
    masksize = mask.size * geo['mask']['resolution'] * 1e-3

    # # Rotation vector (intrinsic-zyx)
    # alpha = geo['mask']['focus']['anglez'] * np.pi / 180
    # beta = geo['mask']['focus']['angley'] * np.pi / 180
    # gamma = geo['mask']['focus']['anglex'] * np.pi / 180
    # rotmat = np.zeros((3, 3), dtype='float32')
    # rotmat[0, 0] = np.cos(alpha) * np.cos(beta)
    # rotmat[0, 1] = np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)
    # rotmat[0, 2] = np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma)
    # rotmat[1, 0] = np.sin(alpha) * np.cos(beta)
    # rotmat[1, 1] = np.sin(alpha) * np.sin(beta) * np.sin(gamma) - np.cos(alpha) * np.cos(gamma)
    # rotmat[1, 2] = np.sin(alpha) * np.sin(beta) * np.cos(gamma) + np.cos(alpha) * np.sin(gamma)
    # rotmat[2, 0] = -np.sin(beta)
    # rotmat[2, 1] = np.cos(beta) * np.sin(gamma)
    # rotmat[2, 2] = np.cos(beta) * np.cos(gamma)

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

    # Mask position
    mc = np.array([
        geo['mask']['focus']['cenx'], 
        geo['mask']['focus']['dist'],
        geo['mask']['focus']['cenz']], dtype='float32') 

    # Scan direction
    rodrot = rotmatrix(geo['scanner']['rot'])
    direction = np.dot(rodrot, geo['scanner']['axis'])

    # Points on the mask for ray tracing
    p1 = np.array([0.5 * masksize - pos, 0, 100], dtype='float32')
    p2 = np.array([0.5 * masksize - pos, 0, -100], dtype='float32')
    p1 = np.dot(rotmat, p1) + mc
    p2 = np.dot(rotmat, p2) + mc
    p1 = p1 + direction * geo['mask']['shift']
    p2 = p2 + direction * geo['mask']['shift']

    # Intersection of the source and the plane 
    intersectionsource1 = intersect2(s1, s2, p0, p1, p2)
    intersectionx = intersectionsource1[0]
    
    # # Magnification
    # _p1 = p1 + direction * geo['scanner']['step']
    # _p2 = p2 + direction * geo['scanner']['step']
    # intersectionsource2 = intersect2(s1, s2, p0, _p1, _p2)
    # magnification = np.sqrt(np.sum(np.power(intersectionsource2 - intersectionsource1, 2)))
    magnification = 1
    
    # Discrete grid along the source beam
    gr = geo['source']['grid']
    grid = np.arange(gr[0], gr[1], gr[2])
    extgrid = np.arange(gr[0] - 1.0, gr[1] + 1.0, gr[2])

    # Stretch factor for magnification
    angpix = np.arctan(p0[0] / p0[1]) 
    angmsk = geo['mask']['focus']['anglez'] * np.pi / 180.
    xxx = np.cos(angpix) / np.cos(angpix + angmsk)
    _sig = stretcharr(sig, magnification * geo['mask']['resolution'] / (1e3 * gr[2] * xxx), order=1)
    # print (factor, xxx, sig.shape, _sig.shape)

    # Index along the beam
    dx = np.argmin(np.abs(extgrid - intersectionx))

    # Sub-pixel adjustment based on intersectionx
    sx = (intersectionx / gr[2] - np.round(intersectionx / gr[2]))
    _sig = ndimage.shift(_sig, sx, order=1)

    # Initialize fine footprint
    fp = np.zeros((len(extgrid), ), dtype='float32')
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
    try:
        a = (sig.size * c - p * d) / (p - sig.size)
        factor =  p * a / (a + c) / sig.size  / (1e3 * gr[2])
    except ZeroDivisionError:
        factor = 1
    return stretcharr(sig, factor * geo['mask']['resolution'], order=1)


def intersect(s1, s2, p0, p1, p2):
    s12 = (s2 - s1) / np.linalg.norm(s2 - s1)
    p01 = (p1 - p0) / np.linalg.norm(p1 - p0) # vector A
    p02 = (p2 - p0) / np.linalg.norm(p2 - p0) # vector B
    vec = np.cross(p01, p02) # Normal vector
    t = np.dot(vec, s1 - p0) / np.dot(-s12, vec)
    intersect = s1[0] + t * s12[0]
    return intersect


def intersect2(s1, s2, p0, p1, p2):
    s12 = (s2 - s1) / np.linalg.norm(s2 - s1)
    p01 = (p1 - p0) / np.linalg.norm(p1 - p0) # vector A
    p02 = (p2 - p0) / np.linalg.norm(p2 - p0) # vector B
    vec = np.cross(p01, p02) # Normal vector
    t = np.dot(vec, s1 - p0) / np.dot(-s12, vec)
    intersect = s1 + t * s12
    return intersect


def pix2pos(index, geo):
    """Returns the coordinate of a given detector pixel index."""
    dx = geo['detector']['shape'][0]
    dy = geo['detector']['shape'][1]    
    sx = geo['detector']['size'][0]
    sy = geo['detector']['size'][1]
    xp = (index[1] - 0.5 * (dx - 1)) * sx / dx
    yp = (index[0] - 0.5 * (dy - 1)) * sy / dy
    zp = 0
    xp += geo['detector']['pos'][0]
    yp += geo['detector']['pos'][1]
    zp += geo['detector']['pos'][2]
    xpoint = np.array([xp, yp, zp])
    rodrot = rotmatrix(geo['detector']['rot'])
    point = np.dot(rodrot, xpoint)
    # return np.array([point[0], point[1], point[2]], dtype='float32')
    return np.array([point[2], point[1], -point[0]], dtype='float32')


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
