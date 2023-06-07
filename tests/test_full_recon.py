import unittest
import h5py
import indices
import dataclasses
import numpy as np
import cold

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
    def _run_cold(cls, config_path, dat, ind):
        file, comp, geo, algo = cold.config(config_path)

        pos, sig, scl, ene = cold.decode(dat, ind, comp, geo, algo)
        dep, lau = cold.resolve(dat, ind, pos, sig, geo, comp)

        test_results = TestResult(
            pos=np.asarray(pos),
            sig=np.asarray(sig),
            scl=np.asarray(scl),
            ene=np.asarray(ene),
            dep=np.asarray(dep),
            lau=np.asarray(lau)
        )
        return test_results


class TestSI(ColdTestbed):
    @classmethod
    def setUpClass(cls):
        cls.test_data, cls.expected = cls._load_test_data('data/test_si_x1800.h5')
        cls.result = cls._run_cold('configs/SiNorcada90_calib_pos2_maskX1800.yml', 
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
        cls.result = cls._run_cold('configs/twin_pristine_1.yml', 
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