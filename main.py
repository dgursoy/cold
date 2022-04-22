#!/usr/bin/env python3

import cold
import fire
import dxchange

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
    file, comp, geo, algo = cold.config(path)
    for m in range(441):
        file['range'] = [m * 201, 201 * (m + 1), 1]
        data, ind = cold.load(file)
        pos, sig, scl = cold.decode(data, ind, comp, geo, algo, debug=debug)
        # dist = cold.calibrate(data, ind, pos, sig, geo, comp)
        # geo['mask']['focus']['dist'] = dist
        dep, lau = cold.resolve(data, ind, pos, sig, geo, comp)

        shape = geo['detector']['shape']
        cold.saveimg('tmp/pos/pos-' + str(m), pos, ind, shape)
        cold.plotarr('tmp/sig/sig-' + str(m), sig, plots=False)
        cold.saveplt('tmp/dep/dep-' + str(m), dep, geo['source']['grid'])
        cold.saveimg('tmp/lau/lau-' + str(m), lau, ind, shape, swap=True)

if __name__ == '__main__':
    fire.Fire(main)