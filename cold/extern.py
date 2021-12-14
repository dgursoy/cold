#!/usr/bin/env python3

import numpy as np
import ctypes


__author__ = "Doga Gursoy"
__copyright__ = "Copyright (c) 2021, UChicago Argonne, LLC."
__docformat__ = 'restructuredtext en'


LIBCOLD = ctypes.cdll.LoadLibrary(
    'build/lib.macosx-10.7-x86_64-3.6/libcold.cpython-36m-darwin.so')



def recpos(data, mask, pos, sig, tau, dec, method):
    """Returns the mask position for a single pixel."""
    if method == 'lsqr':
        LIBCOLD.recposlsqr.restype = ctypes.c_int
        pos = LIBCOLD.recposlsqr(
            as_c_float_p(data),
            as_c_int(data.size),
            as_c_float_p(mask),
            as_c_int(mask.size),
            as_c_int(0), 
            as_c_int(pos), 
            as_c_int(sig), 
            as_c_float(tau), 
            as_c_float(dec))
    elif method == 'maxl':
        pass
    return pos


def recsig(data, mask, pos, sig, tau, dec, method):
    """Returns the signal footprint size for a single pixel."""
    if method == 'lsqr':
        LIBCOLD.recsiglsqr.restype = ctypes.c_int
        sig = LIBCOLD.recsiglsqr(
            as_c_float_p(data),
            as_c_int(data.size),
            as_c_float_p(mask),
            as_c_int(mask.size),
            as_c_int(0), 
            as_c_int(pos), 
            as_c_int(sig), 
            as_c_float(tau), 
            as_c_float(dec))
    elif method == 'maxl':
        pass
    return sig


def rectau(data, mask, pos, sig, tau, dec, method):
    """Returns the signal footprint size for a single pixel."""
    if method == 'lsqr':
        LIBCOLD.rectaulsqr.restype = ctypes.c_float
        tau = LIBCOLD.rectaulsqr(
            as_c_float_p(data),
            as_c_int(data.size),
            as_c_float_p(mask),
            as_c_int(mask.size),
            as_c_int(0), 
            as_c_int(pos), 
            as_c_int(sig), 
            as_c_float(tau), 
            as_c_float(dec))
    elif method == 'maxl':
        pass
    return tau


def recdec(data, mask, pos, sig, tau, dec, method):
    """Returns the signal footprint size for a single pixel."""
    if method == 'lsqr':
        LIBCOLD.recdeclsqr.restype = ctypes.c_float
        tau = LIBCOLD.recdeclsqr(
            as_c_float_p(data),
            as_c_int(data.size),
            as_c_float_p(mask),
            as_c_int(mask.size),
            as_c_int(0), 
            as_c_int(pos), 
            as_c_int(sig), 
            as_c_float(tau), 
            as_c_float(dec))
    elif method == 'maxl':
        pass
    return tau


def footprint(n, alpha, decay):
    """Returns a footprint."""
    sig = np.zeros((n, ), dtype='float32')
    LIBCOLD.footprint.restype = as_c_void_p()
    LIBCOLD.footprint(
        as_c_int(n),
        as_c_float(alpha), 
        as_c_float(decay), 
        as_c_float_p(sig))
    return sig


def convolve(mask, kernel):
    """Returns the convolution of a signal with a gi en kernel."""
    res = np.zeros((mask.size, ), dtype='float32')
    LIBCOLD.convolve.restype = as_c_void_p()
    LIBCOLD.convolve(
        as_c_float_p(mask),
        as_c_int(mask.size),
        as_c_int_p(kernel),
        as_c_int(kernel.size),
        as_c_float_p(res))
    return res

def as_c_float_p(arr):
    c_float_p = ctypes.POINTER(ctypes.c_float)
    return arr.ctypes.data_as(c_float_p)


def as_c_int_p(arr):
    c_int_p = ctypes.POINTER(ctypes.c_int)
    return arr.ctypes.data_as(c_int_p)


def as_c_long_p(arr):
    c_long_p = ctypes.POINTER(ctypes.c_long)
    return arr.ctypes.data_as(c_long_p)


def as_c_int(arr):
    return ctypes.c_int(arr)


def as_c_long(arr):
    return ctypes.c_long(arr)


def as_c_float(arr):
    return ctypes.c_float(arr)


def as_c_void_p():
    return ctypes.POINTER(ctypes.c_void_p)
