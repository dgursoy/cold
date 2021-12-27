#!/usr/bin/env python3

import numpy as np
import logging
from skimage.feature import blob_log
import matplotlib.pyplot as plt
from scipy import signal, ndimage, optimize, linalg



def smooth(img, sigma):
    img = ndimage.gaussian_filter(img, sigma)
    logging.info("Image smoothed.")
    return img


def blobsearch(img, threshold, minsize):
    blobs = blob_log(img, max_sigma=30, num_sigma=10, threshold=5)
    logging.info("Number of blobs found: " + str(blobs.shape[0]))
    return blobs


def sortblobs(blobs, img):
    nblobs = blobs.shape[0]
    arr = np.zeros((nblobs, ), dtype='float32')
    for m in range(nblobs):
        arr[m] = img[int(blobs[m, 0]), int(blobs[m, 1])]
    blobs = blobs[np.argsort(arr)]
    return blobs


def processblobs(blobs):
    peaks = None
    return peaks


def peaksearch():
    # img = smooth(img, 2)
    # blobs = blobsearch(img, 0, 10)
    # blobs = sortblobs(blobs, img)
    blobs = np.load('trunk/blobs.npy')
    img = np.load('trunk/img.npy')
    peaks = processblobs(blobs, img)
    peaks = None
    return peaks