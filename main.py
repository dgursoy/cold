#!/usr/bin/env python3

import cold
import fire


def main(path, debug=False):
    """Runs the reconstruction workflow given parameters 
    in a configuration file.

    Parameters
    ----------
    path: string
        Path of the YAML file with configuration parameters.

    debug: bool
        If True, plots the fitted signals. 

    Returns
    -------
        None
    """

    file, geo, algo = cold.config(path)
    dat, ind = cold.load(file)
    msk, grd = cold.mask(geo['mask'])
    pos, sig, tau = cold.decode(dat, msk, grd, geo, algo)
    dep, lau = cold.resolve(dat, ind, pos, sig, tau, geo)

    shape = geo['detector']['shape']
    cold.saveimg('tmp/pos/pos', pos, ind, shape)
    cold.saveimg('tmp/sig/sig', sig, ind, shape)
    cold.saveimg('tmp/tau/tau', tau, ind, shape)
    cold.saveimg('tmp/lau/lau', lau, ind, shape, swap=True)
    cold.saveplt('tmp/dep/dep', dep, geo['source']['grid'])
    if debug is True:
        cold.plotresults(dat, ind, msk, grd, pos, sig, tau, geo, algo)

if __name__ == '__main__':
    fire.Fire(main)