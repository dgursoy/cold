import unittest
import h5py
import indices
import dataclasses
import numpy as np
import cold
import time

@dataclasses.dataclass
class TestData:
    data: np.ndarray
    ind: np.ndarray

@dataclasses.dataclass
class TestResult:
    pos: np.ndarray
    sig: np.ndarray
    scl: np.ndarray
    ene: np.ndarray
    dep: np.ndarray
    lau: np.ndarray

@dataclasses.dataclass
class ColdConfig:
    file: dict
    comp: dict
    geo: dict
    algo: dict


class ColdTestbed(unittest.TestCase):

    @classmethod
    def _load_test_data(cls, input_fp):
        with h5py.File(input_fp, 'r') as data_f:
            test_data = TestData(
                data=np.asarray(data_f['data']),
                ind=np.asarray(data_f['ind']),
            )

            test_results = TestResult(
                pos=np.asarray(data_f['pos']),
                sig=np.asarray(data_f['sig']),
                scl=np.asarray(data_f['scl']),
                ene=np.asarray(data_f['ene']),
                dep=np.asarray(data_f['dep']),
                lau=np.asarray(data_f['lau'])
            )
        return test_data, test_results

    @classmethod
    def _run_cold(cls, cold_cfg, dat, ind):

        start_time = time.time()
        pos, sig, scl, ene, _ = cold.decode(dat, ind, cold_cfg.comp, cold_cfg.geo, cold_cfg.algo)
        dep, lau = cold.resolve(dat, ind, pos, sig, cold_cfg.geo, cold_cfg.comp)

        test_results = TestResult(
            pos=np.asarray(pos),
            sig=np.asarray(sig),
            scl=np.asarray(scl),
            ene=np.asarray(ene),
            dep=np.asarray(dep),
            lau=np.asarray(lau)
        )
        print(f'Runtime: {time.time() - start_time}')
        return test_results


class TestSI(ColdTestbed):
    @classmethod
    def setUpClass(cls):
        cls.test_data, cls.expected = cls._load_test_data('data/test_si_x1800.h5')
        file, comp, geo, algo = cold.config('configs/SiNorcada90_calib_pos2_maskX1800.yml')
        cold_cfg = ColdConfig(file, comp, geo, algo)
        cls.result = cls._run_cold(cold_cfg, 
                                   cls.test_data.data, 
                                   cls.test_data.ind)

    def test_pos(self):
        self.assertTrue(np.array_equal(self.result.pos, self.expected.pos))

    def test_sig(self):
        self.assertTrue(np.array_equal(self.result.sig, self.expected.sig))

    def test_scl(self):
        self.assertTrue(np.array_equal(self.result.scl, self.expected.scl))

    def test_ene(self):
        self.assertTrue(np.array_equal(self.result.ene, self.expected.ene))

    def test_dep(self):
        self.assertTrue(np.array_equal(self.result.dep, self.expected.dep))

    def test_lau(self):
        self.assertTrue(np.array_equal(self.result.lau, self.expected.lau))


class TestTP(ColdTestbed):
    @classmethod
    def setUpClass(cls):
        cls.test_data, cls.expected = cls._load_test_data('data/test_tp_1.h5')
        file, comp, geo, algo = cold.config('configs/twin_pristine_1.yml')
        cold_cfg = ColdConfig(file, comp, geo, algo)
        cls.result = cls._run_cold(cold_cfg, 
                                   cls.test_data.data, 
                                   cls.test_data.ind)

    def test_pos(self):
        self.assertTrue(np.array_equal(self.result.pos, self.expected.pos))

    def test_sig(self):
        self.assertTrue(np.array_equal(self.result.sig, self.expected.sig))

    def test_scl(self):
        self.assertTrue(np.array_equal(self.result.scl, self.expected.scl))

    def test_ene(self):
        self.assertTrue(np.array_equal(self.result.ene, self.expected.ene))

    def test_dep(self):
        self.assertTrue(np.array_equal(self.result.dep, self.expected.dep))

    def test_lau(self):
        self.assertTrue(np.array_equal(self.result.lau, self.expected.lau))


class TestSIProc(ColdTestbed):
    @classmethod
    def setUpClass(cls):
        cls.test_data, cls.expected = cls._load_test_data('data/test_si_x1800.h5')
        file, comp, geo, algo = cold.config('configs/SiNorcada90_calib_pos2_maskX1800.yml')
        comp['server'] = 'proc'
        cold_cfg = ColdConfig(file, comp, geo, algo)
        cls.result = cls._run_cold(cold_cfg, 
                                   cls.test_data.data, 
                                   cls.test_data.ind)

    def test_pos(self):
        self.assertTrue(np.array_equal(self.result.pos, self.expected.pos))

    def test_sig(self):
        self.assertTrue(np.array_equal(self.result.sig, self.expected.sig))

    def test_scl(self):
        self.assertTrue(np.array_equal(self.result.scl, self.expected.scl))

    def test_ene(self):
        self.assertTrue(np.array_equal(self.result.ene, self.expected.ene))

    def test_dep(self):
        self.assertTrue(np.array_equal(self.result.dep, self.expected.dep))

    def test_lau(self):
        self.assertTrue(np.array_equal(self.result.lau, self.expected.lau))


class TestTPProc(ColdTestbed):
    @classmethod
    def setUpClass(cls):
        cls.test_data, cls.expected = cls._load_test_data('data/test_tp_1.h5')
        file, comp, geo, algo = cold.config('configs/twin_pristine_1.yml')
        comp['server'] = 'proc'
        cold_cfg = ColdConfig(file, comp, geo, algo)
        cls.result = cls._run_cold(cold_cfg, 
                                   cls.test_data.data, 
                                   cls.test_data.ind)

    def test_pos(self):
        self.assertTrue(np.array_equal(self.result.pos, self.expected.pos))

    def test_sig(self):
        self.assertTrue(np.array_equal(self.result.sig, self.expected.sig))

    def test_scl(self):
        self.assertTrue(np.array_equal(self.result.scl, self.expected.scl))

    def test_ene(self):
        self.assertTrue(np.array_equal(self.result.ene, self.expected.ene))

    def test_dep(self):
        self.assertTrue(np.array_equal(self.result.dep, self.expected.dep))

    def test_lau(self):
        self.assertTrue(np.array_equal(self.result.lau, self.expected.lau))


if __name__ == '__main__':
    unittest.main()

if __name__ == '__main__':
    unittest.main()