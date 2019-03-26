import math
import datetime as dt
import re

import numpy as np
import pandas as pd
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import sys

MEW = 3986004.418e8 #m3/s2
Re = 6378.37e3 # m

class PropagationError(Exception):
    def __init__(self, msg):
        self.msg = msg
    pass

class Object:
    """
    Generic class for orbital object
    """
    def __init__(self, tle_str):
        """
        Initializer for object class, based on TLE
        You must supply a tle_str
        :param string tle_str: tle in string format, used to define object params
        """

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

        ## Calculate last perigee time
        timeSincePerigee = self.M/self.n
        self._pt = self.epoch - dt.timedelta(seconds=timeSincePerigee)

        ## Calculate h,P
        self.h = math.sqrt(MEW*self.a*(1-math.pow(self.e,2)))
        self.P = math.pow(self.h,2)/MEW

        ## Initialize sgp4 object
        self.sgpObject = twoline2rv(line1, line2, wgs72)

    def get_teme_state(self, time):
        """
        Return the object teme position at given time
        :param datetime.datetime time: reference time for position
        :return: list pos, list veloc: position in [x,y,z] format [m] and velocity in [vx, vy, vz] [m]
        """
        r_teme, v_teme = self.sgpObject.propagate(
            time.year, time.month, time.day, time.hour, time.minute, time.second
        ) #km, km/s
        r_teme = [x*1000 for x in r_teme] #m
        v_teme = [v*1000 for v in v_teme] #m

        if self.sgpObject.error != 0:
            raise PropagationError(self.sgpObject.error_message)

        return r_teme, v_teme

    def get_eci_state(self, time):
        ## TODO: Actually rotate from TEME to ECI
        r_eci, v_eci = self.get_teme_state(time)

        return r_eci, v_eci

    def generate_trajectory(self, startTime, endTime, steps):
        """
        Generate a trajectory between startTime and endTime in interval time steps
        :param datetime.datetime startTime: start of trajectory (inclusive)
        :param datetime.datetime endTime: end of trajectory (inclusive)
        :param int steps: number of time steps between trajectory points
        :return: define pandas.DataFrame self.trajectory
        """

        deltaTime = int((endTime - startTime).total_seconds())
        step = round(deltaTime/steps)
        times = list(range(0, deltaTime, step))
        times.append(deltaTime)

        self.trajectoryTimes = []
        self.trajectoryPos = np.zeros([len(times), 3])
        self.trajectoryVeloc = np.zeros([len(times), 3])

        for indx, deltaSeconds in enumerate(times):
            time = startTime + dt.timedelta(seconds=deltaSeconds)

            r, v = self.get_eci_state(time)

            self.trajectoryTimes.append(time)
            self.trajectoryPos[indx,:] = r
            self.trajectoryVeloc[indx,:] = v

        combined_rv = np.hstack((self.trajectoryPos, self.trajectoryVeloc))
        self.trajectory = pd.DataFrame(data=combined_rv, columns=['Posx','Posy','Posz','Velx','Vely','Velz'])
        self.trajectory['Times'] = self.trajectoryTimes

    def parse_trajectory(self, indx):
        """
        Given index in the trajectory DataFrame, return r, v in nice form
        :param indx: indx in self._trajectory
        :return: float list r, float list v in format [x,y,z]
        """

        r = [self.trajectory.iloc[indx]['Posx'], self.trajectory.iloc[indx]['Posy'], self.trajectory.iloc[indx]['Posz']]
        v = [self.trajectory.iloc[indx]['Velx'], self.trajectory.iloc[indx]['Vely'], self.trajectory.iloc[indx]['Velz']]

        return r, v

    def is_in_sphere(self,laser_loc, time, sphere_size):
        """
        Given the location of the laser satellite, the time, and sphere size, calculates if obj is in sphere
        :param list laser_loc: location of laser in ECI [x,y,z] [m]
        :param datetime.datetime time: time in UTC
        :param float sphere_size: radius of targeting sphere [m]
        :return: bool in: True if obj is in the sphere
        """

        ## Find object location at time
        closestPoint, vPoint = self.get_eci_state(time)

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