import cold
import h5py
import numpy as np
import indices

def create_results(config_path, output_path, idxs=None):
    file, comp, geo, algo = cold.config(config_path)
    if indices is not None:
        dat, ind = cold.load(file, index=idxs) 
    else:
        dat, ind = cold.load(file, collapsed=True)

    pos, sig, scl, ene = cold.decode(dat, ind, comp, geo, algo)
    dep, lau = cold.resolve(dat, ind, pos, sig, geo, comp)

    print(f'DEBUG: dat {np.sum(dat)}')
    print(f'DEBUG: ind {np.sum(ind)}')
    print(f'DEBUG: pos {np.sum(pos)}')
    print(f'DEBUG: sig {np.sum(sig)}')
    print(f'DEBUG: scl {np.sum(scl)}')
    print(f'DEBUG: ene {np.sum(ene)}')
    print(f'DEBUG: dep {np.sum(dep)}')
    print(f'DEBUG: lau {np.sum(lau)}')

    """
    with h5py.File(output_path, 'w') as data_f:
        data_f.create_dataset('data', data=dat)
        data_f.create_dataset('ind', data=ind)

        data_f.create_dataset('pos', data=pos)
        data_f.create_dataset('sig', data=sig)
        data_f.create_dataset('scl', data=scl)
        data_f.create_dataset('ene', data=ene)
        data_f.create_dataset('dep', data=dep)
        data_f.create_dataset('lau', data=lau)
    """


def build_default_results():
    create_results('configs/twin_pristine_1.yml', 'data/tp_1.h5', indices.TP_INDICIES)
    create_results('configs/SiNorcada90_calib_pos2_maskX1800.yml', 'data/si_x1800.h5', indices.SI_INCICIES)



if __name__ == '__main__':
    build_default_results()