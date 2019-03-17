import unittest
import math
import datetime as dt

import tle_parse

class Test_parse_tle(unittest.TestCase):
    def test_function(self):

        line1 = '1     5U 58002B   19075.71745479 -.00000160  00000-0 -17930-3 0  9998'
        line2 = '2     5  34.2430 201.2113 1845233  88.2396 292.7424 10.84775486155534'
        str = '{}\n{}'.format(line1,line2)

        obj = tle_parse.Object(str)
        sat_num = 5
        epoch = '2019-75 17:13:8.09380'
        epoch = dt.datetime.strptime(epoch,'%Y-%j %H:%M:%S.%f')
        i = math.radians(34.243)
        O = math.radians(201.2113)
        e = 0.1845233
        wp = math.radians(88.2396)
        M = math.radians(292.7424)
        n = 10.84775486*2*math.pi/86400

        self.assertAlmostEqual(sat_num,obj.sat_num,5)
        self.assertEqual(epoch,obj.epoch,5)
        self.assertAlmostEqual(i,obj.i,5)
        self.assertAlmostEqual(O,obj.O,5)
        self.assertAlmostEqual(e,obj.e,5)
        self.assertAlmostEqual(wp,obj.wp,5)
        self.assertAlmostEqual(M,obj.M,5)
        self.assertAlmostEqual(n,obj.n,5)

if __name__ == '__main__':
    unittest.main()
