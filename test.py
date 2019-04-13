## Alex Philpott
## Unit tests for astroUtils

import unittest
import math
import datetime as dt
import platform
import threading

import numpy as np
import pandas as pd
import matlab

from astroUtils import Object, Laser
import astroUtils

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

def assert_vector_almost_equal(testCase,true,test):

    for i in range(0,3):
        testCase.assertAlmostEqual(test[i],true[i],2)

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

    line0 = '0 VANGUARD 1'
    line1 = '1     5U 58002B   19075.71745479 -.00000160  00000-0 -17930-3 0  9998'
    line2 = '2     5  34.2430 201.2113 1845233  88.2396 292.7424 10.84775486155534'
    str = '{}\n{}\n{}'.format(line0, line1, line2)

    timeSincePerigee = M/n
    perigee_time = epoch - dt.timedelta(seconds=timeSincePerigee)

    a = 8620026.476
    h = 57610384007.89
    P = 8326524.503

    def test_function_tle_str(self):
        obj = Object(self.str)
        self.assertAlmostEqual(self.sat_num,obj.satNum,5)
        self.assertEqual(self.epoch,obj.epoch)
        self.assertAlmostEqual(self.i,obj.i,5)
        self.assertAlmostEqual(self.O,obj.O,5)
        self.assertAlmostEqual(self.e,obj.e,5)
        self.assertAlmostEqual(self.wp,obj.wp,5)
        self.assertAlmostEqual(self.M,obj.M,5)
        self.assertAlmostEqual(self.n,obj.n,5)
        self.assertAlmostEqual(self.a,obj.a,3)
        self.assertAlmostEqual(self.h,obj.h,3)
        self.assertAlmostEqual(self.P,obj.P,3)
        self.assertEqual(self.perigee_time, obj._pt)
        self.assertFalse(obj.deb)
        self.assertEqual('VANGUARD 1',obj.satName)

    def test_deb_classification(self):
        line0 = '0 THOR ABLESTAR R/B'
        line1 = '1    47U 60007C   19078.62082020 -.00000047 +00000-0 +14842-4 0  9994'
        line2 = '2    47 066.6649 293.7858 0235398 276.6054 080.8306 14.42023965068951'
        str = '{}\n{}\n{}\n'.format(line0, line1, line2)
        obj = Object(tle_str=str)
        self.assertTrue(obj.deb)

        line0 = '0 ECHO 1 DEB (METAL OBJ)'
        line1 = '1    53U 60009E   19078.45142974 -.00000010 +00000-0 +73023-3 0  9994'
        line2 = '2    53 047.2747 225.1950 0099831 173.3107 186.9037 12.17225727613029'
        str = '{}\n{}\n{}\n'.format(line0, line1, line2)
        obj = Object(tle_str=str)
        self.assertTrue(obj.deb)

        line0 = '0 SOLRAD 1 (GREB)'
        line1 = '1    46U 60007B   19078.71696926 -.00000006 +00000-0 +21930-4 0  9990'
        line2 = '2    46 066.6907 197.8491 0203173 111.6021 250.6868 14.49212027073077'
        str = '{}\n{}\n{}\n'.format(line0, line1, line2)
        obj = Object(tle_str=str)
        self.assertFalse(obj.deb)

class Test_solve_kepler(unittest.TestCase):
    def test_function(self):
        E = math.radians(85)
        e = 0.8
        M = 0.68657

        E_out = astroUtils.solve_kepler(M,e)
        self.assertAlmostEqual(E,E_out,2)

class Test_calc_o(unittest.TestCase):
    def test_function(self):
        o = math.radians(85)
        e = 0.8
        E = 2*math.atan2(math.sqrt((1-e)/(1+e))*math.tan(o/2),1)

        o_out = astroUtils.calc_o(E,e)
        self.assertAlmostEqual(o,o_out,2)

class Test_peri2eci(unittest.TestCase):
    def test_function(self):
        i = math.radians(30)
        O = math.radians(40)
        wp = math.radians(60)

        Q = np.array([[-0.0991,0.89593,0.43301],[-0.94175,-0.22496,0.25],[0.32139,-0.38302,0.86603]])
        Q = np.transpose(Q)
        Q_out = astroUtils.peri2eci(O,i,wp)

        assert_matrix_almost_equal(self,Q,Q_out)

class Test_generate_trajectory(unittest.TestCase):
    def test_parse_trajectory(self):
        line0 = '0 VANGUARD 1'
        line1 = '1     5U 58002B   19075.71745479 -.00000160  00000-0 -17930-3 0  9998'
        line2 = '2     5  34.2430 201.2113 1845233  88.2396 292.7424 10.84775486155534'
        str = '{}\n{}\n{}'.format(line0, line1, line2)

        obj = Object(tle_str=str)
        startTime = dt.datetime.now()
        endTime = startTime + dt.timedelta(days=1)
        steps = 24
        matlabLock = threading.RLock()

        obj.generate_trajectory(startTime, endTime, steps, matlabLock)
        self.assertEqual(obj.trajectory.iloc[0]['Times'], startTime)
        self.assertEqual(obj.trajectory.iloc[-1]['Times'], endTime)
        self.assertEqual(len(obj.trajectory['Times']),25)

        r, v = obj.get_teme_state(startTime)
        r = r.tolist()[0]
        v = v.tolist()[0]
        rout, vout = obj.parse_trajectory(indx=0)
        self.vectors_equal(r,rout)
        self.vectors_equal(v,vout)

        rout, vout = obj.parse_trajectory(time=startTime)
        self.vectors_equal(r, rout)
        self.vectors_equal(v, vout)

        ## Test TEME trajectory
        obj.generate_trajectory(startTime, endTime, steps, matlabLock)
        self.assertEqual(obj.trajectory.iloc[0]['Times'], startTime)
        self.assertEqual(obj.trajectory.iloc[-1]['Times'], endTime)
        self.assertEqual(len(obj.trajectory['Times']), 25)

        r, v = obj.get_teme_state(startTime)
        r = r.tolist()[0]
        v = v.tolist()[0]
        rout, vout = obj.parse_trajectory(indx=0)
        self.vectors_equal(r,rout)
        self.vectors_equal(v,vout)

        rout, vout = obj.parse_trajectory(time=startTime)
        self.vectors_equal(r,rout)
        self.vectors_equal(v,vout)

    def vectors_equal(self,test, true):
        self.assertEqual(len(test), 3)

        self.assertEqual(true[0], test[0])
        self.assertEqual(true[1], test[1])
        self.assertEqual(true[2], test[2])

class Test_update_tle(unittest.TestCase):
    def test_function(self):
        if platform.system() == 'Windows':
            line0 = '0 VANGUARD 1'
            line1 = '1     5U 58002B   19075.71745479 -.00000160  00000-0 -17930-3 0  9998'
            line2 = '2     5  34.2430 201.2113 1845233  88.2396 292.7424 10.84775486155534'
            str = '{}\n{}\n{}'.format(line0, line1, line2)

            obj = Object(tle_str=str)
            matlabLock = threading.RLock()
            obj.satNum = 99999
            r = np.array([-1761336.083,5782317.269,-3609765.529])
            v = np.array([-2316.498,-4430.769,-5686.826])

            tle = obj.update_tle(r, v, matlabLock, dt.datetime.now())
        else:
            print('Skipping Test_update_tle since not on Windows')

class Test_get_JD(unittest.TestCase):
    def test_function(self):
        time = dt.datetime(1996,10,26,14,20,0,0)
        JD = astroUtils.get_jd(time)

        self.assertAlmostEqual(JD,2450383.09722222,8)

class Test_teme2eci(unittest.TestCase):
    def test_function(self):
        tle = '{}\n{}\n{}'.format('TEST','1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753',
                                  '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667')
        obj = Object(tle)
        time = obj.epoch + dt.timedelta(days=3)
        matlabLock = threading.RLock()

        r_teme = [[-9060.47373569], [4658.70952502], [813.68673153]] #km
        v_teme = [[-2.232832783], [-4.110453490], [-3.157345433]] #km/s

        r_eci = [[-9059.9413786],[4659.6972000],[813.9588875]]
        v_eci = [[-2.233348094],[-4.110136162],[-3.157394074]]

        r_eci_out, v_eci_out = astroUtils.teme2eci(np.array(r_teme), np.array(v_teme), time,
                                                   dt.timedelta(seconds=32), matlabLock)

        r_eci = np.array(r_eci).T.tolist()[0]
        v_eci = np.array(v_eci).T.tolist()[0]
        assert_vector_almost_equal(self, r_eci, r_eci_out.T.tolist()[0])
        assert_vector_almost_equal(self, v_eci, v_eci_out.T.tolist()[0])

class Test_rnd(unittest.TestCase):
    def test_function(self):

        rounded = astroUtils.rnd(10,5)
        self.assertEqual(rounded, 10)

        rounded = astroUtils.rnd(12.5,5)
        self.assertEqual(rounded,15)

        rounded = astroUtils.rnd(12,5)
        self.assertEqual(rounded,10)

        rounded = astroUtils.rnd(13,2)
        self.assertEqual(rounded,14)

        rounded = astroUtils.rnd(13,5)
        self.assertEqual(rounded,15)

        rounded = astroUtils.rnd(13,13)
        self.assertEqual(rounded,13)

class Test_get_bounds(unittest.TestCase):
    def test_function(self):

        x = [0,10,20,30,40,50,60,70,80,90,100]
        y = [0,1,2,3,4,5,6,7,8,9,10]
        data = np.zeros([len(y),len(x)])

        xb, yb = astroUtils.get_bounds(data,x,y)
        self.assertEqual(xb[0],0)
        self.assertEqual(xb[1], 100)
        self.assertEqual(yb[0], 0)
        self.assertEqual(yb[1], 10)

        data[3][4] = 5
        xb, yb = astroUtils.get_bounds(data, x, y)
        self.assertEqual(xb[0], 0)
        self.assertEqual(xb[1], 100)
        self.assertEqual(yb[0], 0)
        self.assertEqual(yb[1], 10)

        data[7][8] = 10
        xb, yb = astroUtils.get_bounds(data, x, y)
        self.assertEqual(xb[0], 40)
        self.assertEqual(xb[1], 80)
        self.assertEqual(yb[0], 3)
        self.assertEqual(yb[1], 7)

        data[7][9] = 12
        xb, yb = astroUtils.get_bounds(data, x, y)
        self.assertEqual(xb[0], 40)
        self.assertEqual(xb[1], 90)
        self.assertEqual(yb[0], 3)
        self.assertEqual(yb[1], 7)

        data[10][9] = 12
        xb, yb = astroUtils.get_bounds(data, x, y)
        self.assertEqual(xb[0], 40)
        self.assertEqual(xb[1], 90)
        self.assertEqual(yb[0], 3)
        self.assertEqual(yb[1], 10)

class Test_append(unittest.TestCase):
    def test_function(self):

        foo = pd.DataFrame([],columns=['x','y'])
        bar = pd.DataFrame([[1,2]],columns=['x','y'])

        foo = astroUtils.append(foo,bar)
        self.assertEqual(foo.iloc[0]['x'],1)
        self.assertEqual(foo.iloc[0]['y'],2)
        self.assertEqual(len(foo),1)

        bar = pd.DataFrame([[3, 4]], columns=['x', 'y'])
        foo = astroUtils.append(foo,bar)
        self.assertEqual(foo.iloc[0]['x'], 1)
        self.assertEqual(foo.iloc[0]['y'], 2)
        self.assertEqual(foo.iloc[1]['x'], 3)
        self.assertEqual(foo.iloc[1]['y'], 4)
        self.assertEqual(len(foo), 2)

        bar = pd.Series([5,6], index=['x','y'])
        foo = astroUtils.append(foo, bar)
        self.assertEqual(foo.iloc[0]['x'], 1)
        self.assertEqual(foo.iloc[0]['y'], 2)
        self.assertEqual(foo.iloc[1]['x'], 3)
        self.assertEqual(foo.iloc[1]['y'], 4)
        self.assertEqual(foo.iloc[2]['x'], 5)
        self.assertEqual(foo.iloc[2]['y'], 6)
        self.assertEqual(len(foo), 3)

        foo = pd.DataFrame([], columns=['x', 'y'])
        bar = pd.Series([7,8], index=['x','y'])
        foo = astroUtils.append(foo, bar)
        self.assertEqual(type(foo),pd.DataFrame)
        self.assertEqual(foo.iloc[0]['x'],7)
        self.assertEqual(foo.iloc[0]['y'],8)

class Test_split_array(unittest.TestCase):
    def test_function(self):
        indices = (np.array([444, 449]),)
        split = [np.array([444]), np.array([449])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

        indices = (np.array([444, 449, 450]),)
        split = [np.array([444]), np.array([449, 450])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

        indices = (np.array([1, 2, 5]),)
        split = [np.array([1,2]), np.array([5])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

        indices = (np.array([0, 1, 2, 3, 10]),)
        split = [np.array([0,1,2,3]), np.array([10])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

        indices = (np.array([0, 1, 2, 3, 10, 11]),)
        split = [np.array([0, 1, 2, 3]), np.array([10, 11])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

        indices = (np.array([0, 1, 2, 3, 5, 10, 11]),)
        split = [np.array([0, 1, 2, 3]), np.array([5]), np.array([10, 11])]
        split_out = astroUtils.split_array(indices)
        self.assert_equal(split, split_out)

    def assert_equal(self, true, test):

        self.assertEqual(len(true),len(test))
        for i in range(0, len(true)):
            truei = true[i]
            testi = test[i]

            self.assertEqual(len(truei), len(testi))

            for j in range(0, len(truei)):
                    self.assertEqual(testi[j], truei[j])

class Test_matlab_tle_connection(unittest.TestCase):
    def test_function(self):
        if platform.system() == 'Windows':
            date = [2019,3,27,17,0,0]
            string_out = astroUtils.eng.tle_stk_tester(matlab.double(date),nargout=1)
            expected_str = '{}\n{}'.format('1 23224           19086.70833333  .00003060  00000-0  65926-3 0 00000',
                                           '2 23224 066.7110 273.5175 0238424 271.4862 304.1474 14.42081959000016')

            self.assertEqual(string_out,expected_str)
        else:
            print('Skipping Test_matlab_tle_connection as not running on Windows.')

class Test_is_ready(unittest.TestCase):
    def test_function(self):
        tle = '{}\n{}\n{}'.format('TEST', '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753',
                                  '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667')
        obj = Object(tle)

        laser = Laser(obj, True)
        time = dt.datetime.now()

        self.assertTrue(laser.is_ready(time, 60))

        laser.fire(time,time+dt.timedelta(seconds=60),obj)
        self.assertFalse(laser.is_ready(time,1))
        fireTimes = laser.get_fire_times()
        self.assertEqual(fireTimes.iloc[0]['Start'], time)

        # Single laser fire after
        self.assertFalse(laser.is_ready(time+dt.timedelta(seconds=10),60))
        self.assertFalse(laser.is_ready(time, 60))
        self.assertFalse(laser.is_ready(time + dt.timedelta(seconds=5057), 60))

        self.assertTrue(laser.is_ready(time + dt.timedelta(seconds=16000),60))
        laser.fire(time + dt.timedelta(seconds=16000),time + dt.timedelta(seconds=16060),obj)
        fireTimes = laser.get_fire_times()
        self.assertEqual(fireTimes.iloc[0]['Start'], time)
        self.assertEqual(fireTimes.iloc[1]['Start'], time + dt.timedelta(seconds=16000))

        # Second laser fire after
        self.assertFalse(laser.is_ready(time + dt.timedelta(seconds=16001),1))

        # Laser fire before
        self.assertTrue(laser.is_ready(time - dt.timedelta(seconds=6000), 60))
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=2000), 60))
        self.assertTrue(laser.is_ready(time - dt.timedelta(seconds=2000), 1))
        laser.fire(time - dt.timedelta(seconds=6000), time - dt.timedelta(seconds=5940), obj)
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=2000), 1))
        self.assertFalse(laser.is_ready(time, 1))
        fireTimes = laser.get_fire_times()
        self.assertEqual(fireTimes.iloc[0]['Start'], time - dt.timedelta(seconds=6000))
        self.assertEqual(fireTimes.iloc[1]['Start'], time)
        self.assertEqual(fireTimes.iloc[2]['Start'], time + dt.timedelta(seconds=16000))

        # Negative laser fire tests
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=7000), 60))
        self.assertTrue(laser.is_ready(time - dt.timedelta(seconds=7000), 1))
        self.assertTrue(laser.is_ready(time - dt.timedelta(seconds=12000), 60))

        # Test for fire between the two positive fires
        self.assertTrue(laser.is_ready(time + dt.timedelta(seconds=6000), 60))
        self.assertTrue(laser.is_ready(time + dt.timedelta(seconds=10000), 60))
        self.assertFalse(laser.is_ready(time + dt.timedelta(seconds=10942), 60))
        self.assertFalse(laser.is_ready(time + dt.timedelta(seconds=11060), 60))

        # Second negative laser fire
        self.assertTrue(laser.is_ready(time - dt.timedelta(seconds=14000), 60))
        laser.fire(time - dt.timedelta(seconds=12000), time - dt.timedelta(seconds=12060), obj)
        fireTimes = laser.get_fire_times()
        self.assertEqual(fireTimes.iloc[0]['Start'], time - dt.timedelta(seconds=12000))
        self.assertEqual(fireTimes.iloc[1]['Start'], time - dt.timedelta(seconds=6000))
        self.assertEqual(fireTimes.iloc[2]['Start'], time)
        self.assertEqual(fireTimes.iloc[3]['Start'], time + dt.timedelta(seconds=16000))

        # Negative laser fire time tests
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=8000), 60))
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=6000), 60))
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=6942), 60))
        self.assertFalse(laser.is_ready(time - dt.timedelta(seconds=11060), 60))

        # Check that overlapping fires is false
        laser.fire(time + dt.timedelta(seconds=20000), time + dt.timedelta(seconds=20010), obj)
        self.assertFalse(laser.is_ready(time+dt.timedelta(seconds=20005),1))
        self.assertTrue(laser.is_ready(time + dt.timedelta(seconds=21000), 2))
        laser.fire(time + dt.timedelta(seconds=21000), time + dt.timedelta(seconds=21002), obj)
        self.assertFalse(laser.is_ready(time + dt.timedelta(seconds=21001), 0.5))

        # Check in between times is good if short enough
        self.assertTrue(laser.is_ready(time+dt.timedelta(seconds=6000),60))
        self.assertFalse(laser.is_ready(time+dt.timedelta(seconds=6000),150))

        self.assertTrue(laser.is_ready(time + dt.timedelta(seconds=10000), 60))
        laser.fire(time + dt.timedelta(seconds=10000), time + dt.timedelta(seconds=10060), obj)

        fireTimes = laser.get_fire_times()
        self.assertEqual(fireTimes.iloc[0]['Start'],time - dt.timedelta(seconds=12000))
        self.assertEqual(fireTimes.iloc[1]['Start'],time - dt.timedelta(seconds=6000))
        self.assertEqual(fireTimes.iloc[2]['Start'],time)
        self.assertEqual(fireTimes.iloc[3]['Start'],time+dt.timedelta(seconds=10000))
        self.assertEqual(fireTimes.iloc[4]['Start'],time + dt.timedelta(seconds=16000))
        self.assertEqual(fireTimes.iloc[5]['Start'],time + dt.timedelta(seconds=20000))
        self.assertEqual(fireTimes.iloc[6]['Start'], time + dt.timedelta(seconds=21000))

if __name__ == '__main__':
    unittest.main()
