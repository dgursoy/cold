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
    data, ind = cold.load(file)
    mask = cold.mask(geo['mask'])
    pos, sig = cold.decode(data, mask, geo, algo)
    dep, lau = cold.resolve(data, ind, pos, sig, geo)

    shape = geo['detector']['shape']
    cold.saveimg('tmp/pos/pos', pos, ind, shape)
    cold.plotarr('tmp/sig/sig', sig, plots=False)
    cold.saveplt('tmp/dep/dep', dep, geo['source']['grid'])
    cold.saveimg('tmp/lau/lau', lau, ind, shape, swap=True)
    if debug is True:
        cold.plotresults(data, ind, mask, pos, sig, geo, algo)

if __name__ == '__main__':
    fire.Fire(main)