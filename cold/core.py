#!/usr/bin/env python3


import numpy as np
import logging
from skimage import restoration
from scipy import signal, ndimage, optimize, linalg, interpolate
from scipy.interpolate import BSpline
from cold import pack, saveplt, partition, smooth, collapse, loadsingle, load, saveimg
from skimage.feature import blob_log
import multiprocessing
import os
import warnings
warnings.filterwarnings('ignore')

__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'



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


def ipos(file, mask, geo, algo, threshold=1000):
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
    x = np.arange(0, 2048, 1)
    y = np.arange(0, 2048, 1)
    xx, yy = np.meshgrid(x, y)
    pos0 = interpolate.griddata(ind, pos, (yy, xx), method='nearest')

    # xgrid = np.mgrid[0:2048:2048j, 0:2048:2048j]
    # xflat = xgrid.reshape(2, -1).T
    # RR = interpolate.RBFInterpolator(ind, pos, kernel='thin_plate_spline') # linear thin_plate_spline cubic quintic
    # yflat = RR(xflat)
    # pos0 = yflat.reshape(2048, 2048)

    # Testing
    test = np.zeros(img.shape)
    for m in range(blobs.shape[0]):
        test[points[m, 0], points[m, 1]] = 1
    test = smooth(test, 2)

    # Save
    import dxchange
    dxchange.write_tiff(pos0 +  test / test.max() * 1e3, 'tmp/pos0.tiff')
    dxchange.write_tiff(pos0, 'tmp/pos0.tiff')
    return pos0









def remoterec(data, ind, mask, geo, algo):

    # Number pf chunks
    chunks = len(data)

    # Initialize mask
    mask = invert(mask)

    # Number of total pixels
    npix = 0
    for item in data:
        npix += item.shape[0]

    # Initialize pos
    pos = np.zeros((npix, ), dtype='float32')
    pos = partition(pos, chunks)
    
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
    sig = np.tile(sig, (npix, 1))
    sig = partition(sig, chunks)

    # Initialize scl
    scl = []
    factor = geo['scanner']['step'] / geo['mask']['resolution']
    beta = (90 - geo['scanner']['angle']) * np.pi / 180 
    weighty = 1 / np.cos(geo['mask']['focus']['angley'] * np.pi / 180)
    weightz = np.cos(geo['mask']['focus']['anglez'] * np.pi / 180)
    for n in range(chunks):
        sclscan = np.ones((ind[n].shape[0], ), dtype='float32')
        for m in range(ind[n].shape[0]):
            p0 = pix2pos(ind[n][m], geo) # [<->, dis2det, v^]
            alpha = np.arctan2(p0[0], p0[1]) #+25
            sclscan[m] = np.sin(beta) + np.cos(beta) * np.tan(alpha)
        scl.append(weighty * weightz * factor * geo['scanner']['step'] * sclscan)

    if algo['server'] == 'local':
        # Pack arguments as list and run
        mask = [mask] * chunks
        algo = [algo] * chunks
        args = [data, mask, algo, pos, sig, scl]
        results = runpar(_decode, args, chunks)

        # Unpack results and rescale them
        for m in range(chunks):
            pos[m] = results[m][0]
            sig[m] = results[m][1]

    if algo['server'] == 'remote':
        from funcx.sdk.client import FuncXClient
        import time

        fxc = FuncXClient()
        fxc.version_check()

        ep_info = fxc.get_endpoint_status(algo['functionid'])
        if ep_info['status'] != 'online': 
            raise RuntimeError("End point is not online!")

        fid = fxc.register_function(_remoterec, description=f"Decoder.")

        # Create a batch function tasks
        job = fxc.create_batch()
        for m in range(chunks):
            job.add(
                data[m], mask, algo, pos[m], sig[m], scl[m],
                endpoint_id=algo['functionid'], function_id=fid)

        batch_task_ids = fxc.batch_run(job)
        counter = 0
        while True: 
            batch_task_status = fxc.get_batch_result(batch_task_ids)
            finished_tasks = [s for s in batch_task_status if batch_task_status[s]['status'] == 'success']
            running_tasks = [s for s in batch_task_status if batch_task_status[s]['status'] != 'success']
            if not running_tasks: break
            else: 
                print(f"Total exec. time (sec): {counter}; Tasks are still running: {len(finished_tasks)} / {len(running_tasks)}")
            time.sleep(10) # Wait 2 secs before checking the results
            counter += 10

        # Unpack results and rescale them
        for m, task in zip(range(chunks), finished_tasks): 
            pos[m] = batch_task_status[task]['result'][0]
            sig[m] = batch_task_status[task]['result'][1]
            print(batch_task_status[task]['result'][2])

    return pos, sig


def _remoterec(data, mask, algo, pos, sig, scl):
    from scipy import ndimage
    import cold
    import time
    t = time.time()

    # Init bases
    size = sig.shape[1]
    order = algo['sig']['order']
    scale = algo['sig']['scale']
    grid = np.arange(0, size + 1, dtype='int')
    knots = np.arange(0, size + 1, scale, dtype='int')
    ncoefs = knots.size - order - 1
    bases = np.zeros((size, ncoefs))
    for m in range(ncoefs):
        t = knots[m:m + order + 2] # knots
        b = BSpline.basis_element(t)
        inc = knots[m]
        while knots[m + order + 1] > inc:
            bases[inc, m] = b.integrate(grid[inc], grid[inc + 1])
            inc += 1

    # for m in range(data.shape[0]):
    #     _data = cold.normalize(data[m])
    #     _data = ndimage.zoom(_data, scl[m], order=1)
    #     _pos = cold.posrecon(_data, mask, pos[m], sig[m], algo)
    #     _sig = cold.sigrecon(_data, mask, _pos, sig[m], algo)
    #     pos[m] = _pos
    #     sig[m] = _sig

    elapsedtime = time.time() - t
    return pos, sig, elapsedtime




def pixrecon(data, mask, pos, sig, scale, algo, pos0, lamda, bases):
    import cold
    import time
    from scipy import ndimage
    t = time.time()
    for m in range(data.shape[0]):
        _data = cold.normalize(data[m])
        _data = ndimage.zoom(_data, scale[m], order=1)
        _pos = cold.posrecon(_data, mask, sig[m], algo, pos0[m], lamda)
        _sig = cold.sigrecon(_data, mask, _pos, sig[m], algo, bases)
        pos[m] = _pos
        sig[m] = _sig
    elapsedtime = time.time() - t
    return pos, sig, elapsedtime




def recon(data, ind, mask, geo, algo, pos0=None, lamda=None):

    # Number pf processes
    chunks = len(data)

    # Number of total pixels
    npix = 0
    for item in data:
        npix += item.shape[0]

    scales = initscales(ind, geo, chunks)
    mask = initmask(mask)
    sig = initsig(algo, chunks, npix, 1 / geo['mask']['resolution'])
    pos, pos0, lamda = initpos(chunks, npix, ind, pos0)
    bases = initbases(algo, sig)

    if algo['server'] == 'local':

        # Pack arguments as list and run
        args = packdecodeargs(data, mask, pos, sig, scales, pos0, lamda, algo, bases, chunks)
        results = runpar(_decode, args, chunks)

        # Unpack results and rescale them
        for m in range(chunks):
            pos[m] = results[m][0]
            sig[m] = results[m][1]

    if algo['server'] == 'remote':

        from funcx.sdk.client import FuncXClient
        import time

        fxc = FuncXClient()
        fxc.version_check()

        ep_info = fxc.get_endpoint_status(algo['functionid'])
        if ep_info['status'] != 'online': 
            raise RuntimeError("End point is not online!")

        pixrecon_fid = fxc.register_function(pixrecon, description=f"Decoder.")

        # Create a batch of pi function tasks
        pixrecon_batch = fxc.create_batch()
        for m in range(chunks):
            pixrecon_batch.add(
                data[m], mask, pos[m], sig[m], 
                scales[m], algo, pos0[m], lamda, bases, 
                endpoint_id=algo['functionid'], function_id=pixrecon_fid)

        batch_task_ids = fxc.batch_run(pixrecon_batch)
        counter = 0
        while True: 
            batch_task_status = fxc.get_batch_result(batch_task_ids)
            finished_tasks = [s for s in batch_task_status if batch_task_status[s]['status'] == 'success']
            running_tasks = [s for s in batch_task_status if batch_task_status[s]['status'] != 'success']
            if not running_tasks: break
            else: 
                print(f"Total exec. time (sec): {counter}; Tasks are still running: {len(finished_tasks)} / {len(running_tasks)}")
            time.sleep(10) # Wait 2 secs before checking the results
            counter += 10

        # Unpack results and rescale them
        for m, task in zip(range(chunks), finished_tasks): 
            pos[m] = batch_task_status[task]['result'][0]
            sig[m] = batch_task_status[task]['result'][1]
            print(batch_task_status[task]['result'][2])
         
    return pos, sig, scales, pos0

def initbases(algo, sig):
    bases = baseskernel(sig[0].shape[1], 
                algo['sig']['order'], 
                algo['sig']['scale'])
    return bases

def initscales(ind, geo, chunks):
    args = packscaleargs(ind, geo, chunks)
    scales = runpar(_getscales, args, chunks)
    return scales



def decode(data, ind, mask, geo, algo, pos0=None):
    """Decodes the position and pixel footprints and their positions
    on the mask using coded measurement data."""

    # Number pf processes
    chunks = len(data)

    # Number of total pixels
    npix = 0
    for item in data:
        npix += item.shape[0]

    # Initialize scalings
    args = packscaleargs(ind, geo, chunks)
    scales = runpar(_getscales, args, chunks)

    # Initialize parameters
    mask = initmask(mask)
    sig = initsig(algo, chunks, npix, 1 / geo['mask']['resolution'])

    # Initial position
    pos = initpos(chunks, npix)
    if pos0 is None:
        pos0 = pos
        lamda = 0
    else:
        pos0 = collapse(pos0, pack(ind))
        pos0 = partition(pos0, 8)
        lamda = 1e-5

    # Pack arguments as list and run
    args = packdecodeargs(data, mask, pos, sig, scales, pos0, lamda, algo, chunks)
    results = runpar(_decode, args, chunks)

    # Unpack results and rescale them
    for m in range(chunks):
        pos[m] = results[m][0]
        sig[m] = results[m][1]
    return pos, sig, scales, pos0


def packscaleargs(ind, geo, chunks):
    geo = [geo] * chunks
    return [ind, geo]


def unpackscaleargs(args):
    ind = args[0]
    geo = args[1]
    return ind, geo


def _getscales(args):
    ind, geo = unpackscaleargs(args)
    # Scaling factor
    factor = stretchfactor(
        geo['mask']['resolution'], 
        geo['scanner']['step'])
        
    beta = (90 - geo['scanner']['angle']) * np.pi / 180 
    scalesscan = np.ones((ind.shape[0], ), dtype='float32')
    for m in range(ind.shape[0]):
        p0 = pix2pos(ind[m], geo) # [<->, dis2det, v^]
        alpha = np.arctan2(p0[0], p0[1]) #+25
        scalesscan[m] = np.sin(beta) + np.cos(beta) * np.tan(alpha)

    weight = 1 / np.cos(geo['mask']['focus']['angley'] * np.pi / 180)
    scalesy = weight * np.ones((ind.shape[0], ), dtype='float32')

    weight = np.cos(geo['mask']['focus']['anglez'] * np.pi / 180)
    scalesz = weight * np.ones((ind.shape[0], ), dtype='float32')
    return scalesy * scalesz * scalesscan * factor * geo['scanner']['step']


def stretchlist(data, scales):
    data = data.copy()
    chunks = len(data)
    for m in range(chunks):
        data[m] = stretcharr(data[m], scales)
    return data


def stretcharr(arr, factor, order=1):
    arr = arr.copy()
    if len(arr.shape) == 2:
        zoom = [1, factor]
    elif len(arr.shape) == 1:
        zoom = factor
    arr = ndimage.zoom(arr, zoom, order=order)
    return arr / factor
        
    
def initdata(data, scales):
    return stretchlist(data, scales)


def initmask(mask):
    return invert(mask)


def initpos(chunks, npix, ind, pos0=None):
    pos = np.zeros((npix, ), dtype='float32')
    pos = partition(pos, chunks)
    if pos0 is None:
        pos0 = pos
        lamda = 0
    else:
        pos0 = collapse(pos0, pack(ind))
        pos0 = partition(pos0, 8)
    return pos, pos0, lamda


def initsig(algo, chunks, npix, factor):
    maxsize = algo['sig']['init']['maxsize']
    avgsize = algo['sig']['init']['avgsize']
    sig = np.zeros(maxsize, dtype='float32')
    first = int((maxsize - 1) * 0.5 - avgsize * 0.5)
    sig[first:first+avgsize] = footprnt(avgsize)
    sig = stretcharr(sig, factor)
    sig = np.tile(sig, (npix, 1))
    return partition(sig, chunks)


def packdecodeargs(data, mask, pos, sig, scales, pos0, lamda, algo, bases, chunks):
    mask = [mask] * chunks
    algo = [algo] * chunks
    lamda = [lamda] * chunks
    bases = [bases] * chunks
    return [data, mask, pos, sig, scales, pos0, lamda, algo, bases]


def unpackdecodeargs(args):
    data = args[0]
    mask = args[1]
    pos = args[2]
    sig = args[3]
    scales = args[4]
    pos0 = args[5]
    lamda = args[6]
    algo = args[7]
    bases = args[8]
    return data, mask, pos, sig, scales, pos0, lamda, algo, bases


def stretchfactor(resolution, step):
    factor = step / resolution
    return factor


def _decode(args):
    data, mask, pos, sig, scales, pos0, lamda, algo, bases = unpackdecodeargs(args)
    npix = data.shape[0]

    for m in range(npix):
        pos[m], sig[m] = pixdecode(
            data[m], mask, pos[m], sig[m], scales[m], pos0[m], lamda, algo, bases)
        logging.info('Pixel decoded: ' +
            str(m) + '/' + str(npix - 1) + 
            ' pos=' + str(pos[m].squeeze()) + 
            ' scales=' + str(scales[m].squeeze()))
    return pos, sig


def pixdecode(data, mask, pos, sig, scale, pos0, lamda, algo, bases):
    """The main function for decoding pixel data."""
    data = normalize(data)
    data = ndimage.zoom(data, scale, order=1)
    pos = posrecon(data, mask, sig, algo, pos0, lamda)
    sig = sigrecon(data, mask, pos, sig, algo, bases)
    return pos, sig


# def pixrecon(data, mask, pos, sig, scale, algo, pos0, lamda, bases):
#     from scipy import signal, ndimage, optimize, linalg
#     import numpy as np
#     import cold

#     for p in range(data.shape[0]):
#         _data = data[p]
#         _data = cold.normalize(_data)
#         _data = ndimage.zoom(_data, scale[p], order=1)
#         _sig = sig[p]

#         # position recon
#         sim = signal.convolve(mask, _sig, 'same')
#         costsize = sim.size - _data.size
#         cost = np.zeros((costsize), dtype='float32')
#         for m in range(costsize):
#             if algo['pos']['method'] == 'lsqr':
#                 cost[m] = np.sum(np.power(sim[m:m+_data.size] - _data, 2)) + lamda * np.sum(np.power(m - pos0[p], 2))
#         try:
#             _pos = np.where(cost.min() == cost)[0][0]
#         except IndexError:
#             pass

#         # signal recon
#         first = int((_sig.size - 1) / 2)
#         last = int(mask.size + first - _sig.size - _data.size)
#         if _pos > first and _pos < last:
#             kernel = np.zeros((_data.size, _sig.size), dtype='float32')
#             for m in range(_data.size):
#                 begin = _pos - first + m - 1
#                 end =  begin + _sig.size
#                 kernel[m] = mask[begin:end]
#             if algo['sig']['method'] == 'splines':
#                 coefs = optimize.nnls(np.dot(kernel, bases), _data)[0][::-1]
#                 _sig = np.dot(bases, coefs)
#             if algo['sig']['method'] == 'nnls':
#                 _sig = optimize.nnls(kernel, _data)[0][::-1]
#             if algo['sig']['method'] == 'pinv':
#                 ikernel = linalg.pinv(kernel, algo['sig']['init']['atol'])
#                 _sig = np.dot(ikernel, _data)[::-1]
#             _sig /= _sig.sum()
#         else:
#             _sig *= 0

#         pos[p] = _pos
#         sig[p] = _sig

#     return pos, sig



def pixrecon(data, mask, pos, sig, scale, algo, pos0, lamda, bases):
    import cold
    import time
    from scipy import ndimage
    t = time.time()
    for m in range(data.shape[0]):
        _data = cold.normalize(data[m])
        _data = ndimage.zoom(_data, scale[m], order=1)
        _pos = cold.posrecon(_data, mask, sig[m], algo, pos0[m], lamda)
        _sig = cold.sigrecon(_data, mask, _pos, sig[m], algo, bases)
        pos[m] = _pos
        sig[m] = _sig
    elapsedtime = time.time() - t
    return pos, sig, elapsedtime


def posrecon(data, mask, pos, sig, algo):
    sim = signal.convolve(mask, sig, 'same')
    costsize = sim.size - data.size
    cost = np.zeros((costsize), dtype='float32')
    for m in range(costsize):
        if algo['pos']['method'] == 'lsqr':
            cost[m] = (np.sum(np.power(sim[m:m+data.size] - data, 2)) + 
                algo['pos']['regpar'] * np.sum(np.power(m - pos, 2)))
    try:
        _pos = np.where(cost.min() == cost)[0][0]
    except IndexError:
        pass
    return _pos


def sigrecon(data, mask, pos, sig, algo, bases):
    first = int((sig.size - 1) / 2)
    last = int(mask.size + first - sig.size - data.size)
    if pos > first and pos < last:
        kernel = np.zeros((data.size, sig.size), dtype='float32')
        for m in range(data.size):
            begin = pos - first + m - 1
            end =  begin + sig.size
            kernel[m] = mask[begin:end]
        coefs = optimize.nnls(np.dot(kernel, bases), data)[0][::-1]
        sig = np.dot(bases, coefs)
        sig /= sig.sum()
    else:
        sig *= 0
    return sig


def baseskernel(size, order, scale):
    # Original grid [must be integers]
    grid = np.arange(0, size + 1, dtype='int')

    # Initial knots [must be integers]
    knots = np.arange(0, size + 1, scale, dtype='int')

    # Number of B-splines
    ncoefs = knots.size - order - 1

    bases = np.zeros((size, ncoefs))
    for m in range(ncoefs):

        # B-splines define
        t = knots[m:m + order + 2] # knots
        b = BSpline.basis_element(t)

        inc = knots[m]
        while knots[m + order + 1] > inc:
            bases[inc, m] = b.integrate(grid[inc], grid[inc + 1])
            inc += 1
    return bases


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


def plotresults(dat, ind, msk, pos, sig, scales, geo, algo):
    _dat = pack(dat)
    _ind = pack(ind)
    _pos = pack(pos)
    _sig = pack(sig)
    _scales = pack(scales)
    _msk = msk.copy()
    npix = _ind.shape[0]
    for m in range(npix):
        if _pos[m] > 0:
            if algo['pos']['method'] == 'lsqr':
                plotlsqr(
                    _dat[m], _ind[m], _msk, 
                    _pos[m], _sig[m], _scales[m], m)
            logging.info("Saved: tmp/plot/plot-" + str(m) + ".png")


def plotlsqr(dat, ind, msk, pos, sig, scale, id):
    import matplotlib.pyplot as plt 
    from scipy import ndimage   
    # dat = ndimage.zoom(dat, scale)

    plt.figure(figsize=(8, 4))

    # dat = dat[0:120]
    dat = normalize(dat)
    dat = ndimage.zoom(dat, scale, order=1)
    _sig = sig / sig.max()

    msk = invert(msk)
    _dat = np.zeros(msk.shape)
    _dat[int(pos):int(pos+dat.size)] = dat

    plt.subplot(311)
    plt.title(str(id) + ' ' + ' ind=' + str(ind) + ' pos=' + str(pos))
    # plt.plot(msk[300:6060])
    # plt.plot(_dat[300:6060], 'r')
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

    # plt.figure(figsize=(3, 1))
    # dat = dat[0:120]
    # dat = normalize(dat)
    # plt.plot(dat, 'r')
    # plt.grid('on')
    # plt.ylim((-0.5, 1.5))

    plt.tight_layout()
    if not os.path.exists('tmp/plot'):
        os.makedirs('tmp/plot')
    plt.savefig('tmp/plot/plot-' + str(id) + '.png', dpi=240)
    plt.close()


def calibrate(data, ind, pos, sig, shift, geo):
    # Number pf processes
    chunks = len(data)

    # Initialize 
    dist = np.arange(*geo['mask']['calibrate']['dist'])

    maximum = 0
    optdist = 0
    for k in range(len(dist)):
        # Pack arguments as list and run
        geo['mask']['focus']['dist'] = dist[k]
        args = packcalibrateargs(data, ind, pos, sig, shift, geo, chunks)
        depth = runpar(_calibrate, args, chunks)

        # Save results
        saveplt('tmp/dep-dist/dep-' + str(k), depth, geo['source']['grid'])
        if np.sqrt(np.sum(np.power(depth, 2))) > maximum:
            maximum = np.sqrt(np.sum(np.power(depth, 2)))
            optdist = dist[k]
    logging.info("Optimum distance: " + str(optdist))
    return optdist


def packcalibrateargs(data, ind, pos, sig, shift, geo, chunks):
    geo = [geo] * chunks
    shift = [shift] * chunks
    return [data, ind, pos, sig, shift, geo]


def unpackcalibrateargs(args):
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    shift = args[4]
    geo = args[5]
    return data, ind, pos, sig, shift, geo


def _calibrate(args):
    # Unpack arguments
    data, ind, pos, sig, shift, geo = unpackcalibrateargs(args)

    # Source beam
    grid = np.arange(*geo['source']['grid'])

    # Initializations
    depth = np.zeros((len(grid), ), dtype='float32')

    # Number of pixels to be processed
    npix = data.shape[0]

    for m in range(npix):
        # Calculate the signal footprint along the beam
        fp = _footprint(data[m], ind[m], pos[m], sig[m], shift, geo)
        depth += fp
    return depth


def resolve(data, ind, pos, sig, geo, shift=0):
    """Resolves depth information."""
    # Number pf processes
    chunks = len(data)

    # Pack arguments as list and run
    args = packresolveargs(data, ind, pos, sig, shift, geo, chunks)
    results = runpar(_resolve, args, chunks)

    # Unpack results
    depth = [None] * chunks
    laue = [None] * chunks
    for m in range(chunks):
        depth[m] = results[m][0]
        laue[m] = results[m][1]
    return depth, laue


def packresolveargs(data, ind, pos, sig, shift, geo, chunks):
    shift = [shift] * chunks
    geo = [geo] * chunks
    return [data, ind, pos, sig, shift, geo]


def unpackresolveargs(args):
    data = args[0]
    ind = args[1]
    pos = args[2]
    sig = args[3]
    shift = args[4]
    geo = args[5]
    return data, ind, pos, sig, shift, geo


def _resolve(args):
    # Unpack arguments
    data, ind, pos, sig, shift, geo = unpackresolveargs(args)
    
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
    # p1 = np.array([-pos * 1e-3 + xx, yy, 100 + zz], dtype='float32')
    # p2 = np.array([-pos * 1e-3 + xx, yy, -100 + zz], dtype='float32')
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
    # p1[1], p1[2] = rotpoint2d(p1[1], p1[2], anglex)
    # p2[1], p2[2] = rotpoint2d(p2[1], p2[2], anglex)
    # p1[0], p1[2] = rotpoint2d(p1[0], p1[2], angley)
    # p2[0], p2[2] = rotpoint2d(p2[0], p2[2], angley)
    # p1[0], p1[1] = rotpoint2d(p1[0], p1[1], anglez)
    # p2[0], p2[1] = rotpoint2d(p2[0], p2[1], anglez)
    # p1[0] -= shift
    # p2[0] -= shift

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
