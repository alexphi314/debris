import math
import datetime as dt
import re

import numpy as np

MEW = 3986004.418e8 #m3/s2
Re = 6378.37e3 # m

class Object:
    """
    Generic class for orbital object
    """
    def __init__(self, tle_str = None, kepElems = None):
        """
        Initializer for object class, based on TLE
        You can supply a tle_str or a dict of Keplerian elements. If both are given, the tle string is used
        :param string tle_str: tle in string format, used to define object params
        :param dict kepElems: dictionary with entries {'i','O','e','wp','o','a','epoch','satNum','deb', 'satName'}
                              corresponding to the keplerian orbit elements {Inclin [rad], RAAN [rad], eccent,
                              arg of perigee [rad], true anomaly [rad], semi-major axis [m]} and epoch of the object,
                              satellite number, if it is debris [bool], and satellite name
        """

        if kepElems is not None and tle_str is None:
            self.i = kepElems['i']
            self.O = kepElems['O']
            self.e = kepElems['e']
            self.wp = kepElems['wp']
            o = kepElems['o']
            self.a = kepElems['a']
            self.n = math.sqrt(MEW/math.pow(self.a,3))
            E = 2*math.atan(math.sqrt((1-self.e)/(1+self.e))*math.tan(o/2))
            self.M = E - self.e*math.sin(E)

            if self.M > 2*math.pi:
                self.M = self.M % 2*math.pi

            while self.M < 0:
                self.M += 2*math.pi

            self.epoch = kepElems['epoch']
            self.satNum = kepElems['satNum']
            self.deb = kepElems['deb']
            self.satName = kepElems['satName']
        else:
            lines = tle_str.split('\n')
            line0 = lines[0]
            line1 = lines[1]
            line2 = lines[2]

            self.satName = line0[2:].strip()
            rb_re = re.compile('\S*R/B\S*',re.IGNORECASE)
            deb_re = re.compile('\S*DEB\S*',re.IGNORECASE)

            if rb_re.search(self.satName) or deb_re.search(self.satName):
                self.deb = True
            else:
                self.deb = False

            self.satNum = int(line1[2:7])
            epoch_year = line1[18:20]
            epoch_day = line1[20:32]
            self.i = math.radians(float(line2[8:16]))
            self.O = math.radians(float(line2[17:25]))
            self.e = float("." + line2[26:33])
            self.wp = math.radians(float(line2[34:42]))
            self.M = math.radians(float(line2[43:51]))
            self.n = float(line2[52:63]) * 2 * math.pi / 86400  # rad/s
            self.a = math.pow(MEW / math.pow(self.n, 2), float(1 / 3))  # m

            ## Calculate TLE epoch time
            if int(epoch_year) > int(dt.datetime.now().strftime('%y')):
                year = "19" + str(epoch_year)
            else:
                year = "20" + str(epoch_year)

            frac, doy = math.modf(float(epoch_day))
            frac, hour = math.modf(frac * 24)
            frac, min = math.modf(frac * 60)
            frac, sec = math.modf(frac * 60)

            if doy < 10:
                doy = "00" + str(int(doy))
            elif doy < 100:
                doy = "0" + str(int(doy))
            else:
                doy = str(int(doy))

            epoch = '{}-{} {}:{}:{}.{}'.format(year, doy, int(hour), int(min), int(sec), str(frac)[2:6])
            self.epoch = dt.datetime.strptime(epoch, '%Y-%j %H:%M:%S.%f')

        ## Calculate perigee time
        timeSincePerigee = self.M/self.n
        self._pt = self.epoch - dt.timedelta(seconds=timeSincePerigee)

        ## Calculate h,P
        self.h = math.sqrt(MEW*self.a*(1-math.pow(self.e,2)))
        self.P = math.pow(self.h,2)/MEW

        ## Define peri2eci rotation matrix
        ##TODO: When updating raan and wp, this must be calculated in loop
        self._Q_eci_p = peri2eci(self.O, self.i, self.wp)

        self._intervalTime = 3*60
        self.define_trajectory()

    def define_trajectory(self):
        """
        Create a trajectory vector [[x,y,z],...] in ECI of satellite position throughout 1 orbit
        :return: Define self.trajectory
        """

        self._T = 2*math.pi/self.n
        endTime = int(self._T//self._intervalTime*self._intervalTime) #round T to nearest interval

        times = range(0, endTime+self._intervalTime, self._intervalTime)
        self._trajectory = np.zeros([len(times),3])
        for i in range(0,len(times)):
            ## TODO: Include change in RAAN and wp due to J2 perturbations
            ## Update M, E, and o (true anomaly)
            time = times[i]
            timeSincePerigee = time
            M = self.n*timeSincePerigee
            E = solve_kepler(M, self.e)
            o = calc_o(E, self.e)

            ## Calculate r_eci
            r = self.P / (1 + self.e*math.cos(o))
            r_peri = np.array([[r*math.cos(o)],[r*math.sin(o)],[0]])
            r_eci = np.matmul(self._Q_eci_p, r_peri)

            self._trajectory[i,:] = np.transpose(r_eci)

    def get_eci_pos(self, time):
        """
        Return the object ECI position at given time
        :param datetime time: reference time for position
        :return: list pos: position in [x,y,z] format [m]
        """
        allTimeSincePerigee = (time - self._pt).total_seconds()  # absolute time since last perigee passage

        pt = self._pt - dt.timedelta(seconds=self._T)
        while allTimeSincePerigee < 0:
            allTimeSincePerigee = (time - pt).total_seconds()
            pt -= dt.timedelta(seconds=self._T)

        orbitTimeSincePerigee = allTimeSincePerigee % self._T  # time within 1 orbit since perigee passage
        closestTimeIndex = int(orbitTimeSincePerigee // self._intervalTime)
        pos = self._trajectory[closestTimeIndex, :]

        return pos

    def is_in_sphere(self,laser_loc, time, sphere_size):
        """
        Given the location of the laser satellite, the time, and sphere size, calculates if obj is in sphere
        :param list laser_loc: location of laser in ECI [x,y,z] [m]
        :param datetime time: time in UTC
        :param float sphere_size: radius of targeting sphere [m]
        :return: bool in: True if obj is in the sphere
        """

        ## Find object location at time
        closestPoint = self.get_eci_pos(time)

        ## Compute Distance
        dx = laser_loc[0] - closestPoint[0]
        dy = laser_loc[1] - closestPoint[1]
        dz = laser_loc[2] - closestPoint[2]
        d = math.sqrt(math.pow(dx,2) + math.pow(dy,2) + math.pow(dz,2))

        if d < sphere_size:
            return True

        return False

def solve_kepler(M, e):
    """
    Solve Kepler's Equation for E
    :param float M: mean anomaly [rad]
    :param float e: eccentricity
    :return: float E: eccentric anomaly [rad]
    """

    ## Define equations
    f = lambda x: x - e * math.sin(x) - M
    f_prime = lambda x: 1 - e * math.cos(x)

    ## Pick guess
    if M < math.pi:
        E = M + e / 2
    else:
        E = M - e / 2

    ## Loop until we are close to the root
    ratio = f(E) / f_prime(E)
    while abs(ratio) > 1e-8:
        E -= ratio
        ratio = f(E) / f_prime(E)

    if E > 2.0 * math.pi:
        two_pi = 2.0 * math.pi
        rem = E % two_pi
        E = rem

    while E < 0:
        E += 2.0 * math.pi

    return E

def calc_o(E, e):
    """
    Given eccentric anomaly and eccentricity, calculate true anomaly

    :param float E: eccentric anomaly [rad]
    :param float e: eccentricity
    :return: float o: true anomaly [rad]
    """

    e_coef = math.sqrt((1 - e) / (1 + e))
    o = 2.0 * math.atan(math.tan(E / 2.0) / e_coef)

    if o < 0:
        o += 2.0 * math.pi

    return o

def peri2eci(O, i, wp):
    """
    Given RAAN, inclination, argument of perigee, return O_eci_p
    :param float O: RAAN [rad]
    :param float i: inclination [rad]
    :param float wp: argument of perigee [rad]
    :return: array O_eci_p: DCM from perifocal frame to ECI frame
    """

    q11 = -math.sin(O) * math.cos(i) * math.sin(wp) + math.cos(O) * math.cos(wp)
    q12 = -math.sin(O) * math.cos(i) * math.cos(wp) - math.cos(O) * math.sin(wp)
    q13 = math.sin(O) * math.sin(i)
    q21 = math.cos(O) * math.cos(i) * math.sin(wp) + math.sin(O) * math.cos(wp)
    q22 = math.cos(O) * math.cos(i) * math.cos(wp) - math.sin(O) * math.sin(wp)
    q23 = -math.cos(O) * math.sin(i)
    q31 = math.sin(i) * math.sin(wp)
    q32 = math.sin(i) * math.cos(wp)
    q33 = math.cos(i)

    Q = np.array([[q11, q12, q13], [q21, q22, q23], [q31, q32, q33]])
    return Q

def parse_catalog(format):
    """
    Open the txt file with all TLEs in the sat cat, create an Object for each tle
    :param int format: 2 or 3 corresponding to 2 line or 3 line format
    :return: objects list, each index is an Object in the sat cat
    """
    if format == 2:
        file = 'Data/satcat.txt'
    else:
        file = 'Data/satcat_3le.txt'

    with open(file,'r') as f:
        lines = f.readlines()
        objects = [None]*int(len(lines)/format)
        for i in range(0, len(lines), format):
            tle = lines[i:i+format]

            if format == 2:
                str = '{}{}'.format(tle[0],tle[1])
            else:
                str = '{}{}{}'.format(tle[0], tle[1], tle[2])

            objects[int(i/format)] = Object(tle_str=str)

    return objects