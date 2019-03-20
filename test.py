import unittest
import math
import datetime as dt

import numpy as np

import tle_parse

def assert_matrix_almost_equal(testCase,true,test):
    ## Purpose: Assert that all elements of input matrices are equal
    ##
    ## Inputs:
    ##   True array
    ##   Test array

    r,c = true.shape
    r_test,c_test = test.shape

    testCase.assertEqual(r,r_test)
    testCase.assertEqual(c,c_test)

    for i in range(0,r):
        for j in range(0,c):
            testCase.assertAlmostEqual(test[i][j],true[i][j],2)

class Test_object_init(unittest.TestCase):
    sat_num = 5
    epoch = '2019-75 17:13:8.09380'
    epoch = dt.datetime.strptime(epoch, '%Y-%j %H:%M:%S.%f')
    i = math.radians(34.243)
    O = math.radians(201.2113)
    e = 0.1845233
    wp = math.radians(88.2396)
    M = math.radians(292.7424)
    n = 10.84775486 * 2 * math.pi / 86400

    line1 = '1     5U 58002B   19075.71745479 -.00000160  00000-0 -17930-3 0  9998'
    line2 = '2     5  34.2430 201.2113 1845233  88.2396 292.7424 10.84775486155534'
    str = '{}\n{}'.format(line1, line2)

    def test_function_tle_str(self):
        obj = tle_parse.Object(self.str)
        self.assertAlmostEqual(self.sat_num,obj.satNum,5)
        self.assertEqual(self.epoch,obj.epoch)
        self.assertAlmostEqual(self.i,obj.i,5)
        self.assertAlmostEqual(self.O,obj.O,5)
        self.assertAlmostEqual(self.e,obj.e,5)
        self.assertAlmostEqual(self.wp,obj.wp,5)
        self.assertAlmostEqual(self.M,obj.M,5)
        self.assertAlmostEqual(self.n,obj.n,5)

    def test_function_kepElems(self):
        E = tle_parse.solve_kepler(self.M, self.e)
        o = tle_parse.calc_o(E,self.e)
        a = math.pow(tle_parse.MEW,1/3)/math.pow(self.n,2/3)
        kepElems = {
            'i': self.i,
            'O': self.O,
            'e': self.e,
            'wp': self.wp,
            'epoch': self.epoch,
            'o': o,
            'a': a,
            'satNum': self.sat_num,
            'deb': False
        }
        obj = tle_parse.Object(kepElems=kepElems)
        tle_obj = tle_parse.Object(tle_str=self.str)

        self.assertEqual(obj.satNum, tle_obj.satNum)
        self.assertEqual(obj.epoch, tle_obj.epoch)
        self.assertEqual(obj.i, tle_obj.i)
        self.assertEqual(obj.O, tle_obj.O)
        self.assertEqual(obj.e, tle_obj.e)
        self.assertEqual(obj.wp, tle_obj.wp)
        self.assertAlmostEqual(obj.M, tle_obj.M, 5)
        self.assertAlmostEqual(obj.n, tle_obj.n, 5)
        self.assertEqual(obj.a, tle_obj.a)
        self.assertFalse(obj.deb)

class Test_solve_kepler(unittest.TestCase):
    def test_function(self):
        E = math.radians(85)
        e = 0.8
        M = 0.68657

        E_out = tle_parse.solve_kepler(M,e)
        self.assertAlmostEqual(E,E_out,2)

class Test_calc_o(unittest.TestCase):
    def test_function(self):
        o = math.radians(85)
        e = 0.8
        E = 2*math.atan2(math.sqrt((1-e)/(1+e))*math.tan(o/2),1)

        o_out = tle_parse.calc_o(E,e)
        self.assertAlmostEqual(o,o_out,2)

class Test_peri2geo(unittest.TestCase):
    def test_function(self):
        i = math.radians(30)
        O = math.radians(40)
        wp = math.radians(60)

        Q = np.array([[-0.0991,0.89593,0.43301],[-0.94175,-0.22496,0.25],[0.32139,-0.38302,0.86603]])
        Q = np.transpose(Q)
        Q_out = tle_parse.peri2eci(O,i,wp)

        assert_matrix_almost_equal(self,Q,Q_out)

if __name__ == '__main__':
    unittest.main()
