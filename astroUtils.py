## Alex Philpott
## Utility file defining Object and Laser classes. Also defines various functions used throughout the debSim file.
## Separated to enable easy unit testing.

import math
import datetime as dt
import re
import os

import numpy as np
import pandas as pd
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
import matlab.engine

MEW = 3986004.418e8 #m3/s2
Re = 6378.37e3 # m

## 4/2/19
xp = 48.7e-3 #arcsec
yp = 384.1e-3 #arcsec

xp = math.radians(xp*(0.000028/0.1))
yp = math.radians(yp*(0.000028/0.1))
delta_ut1 = dt.timedelta(seconds=-120.16e-3)

delta_at = dt.timedelta(seconds=37)

## Import MATLAB
eng = matlab.engine.start_matlab()
eng.addpath(os.getcwd()+'/MATLAB')

class PropagationError(Exception):
    def __init__(self, msg, num, time):
        """
        Init error

        :param string msg: error message
        :param int num: error number
        :param datetime.datetime time: trajectory time where error occurred
        """
        self.msg = msg
        self.num = num
        self.time = time
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

        self.a_s = []
        self.i_s = []

        self.tle = tle_str
        self.parse_tle()

    def parse_tle(self):
        """
        Parse the input tle and update variables

        :return:
        """

        lines = self.tle.split('\n')
        line0 = lines[0]
        line1 = lines[1]
        line2 = lines[2]

        self.satName = line0[2:].strip()
        rb_re = re.compile('\S*R/B\S*', re.IGNORECASE)
        deb_re = re.compile('\S*DEB\S*', re.IGNORECASE)

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
        timeSincePerigee = self.M / self.n
        self._pt = self.epoch - dt.timedelta(seconds=timeSincePerigee)

        ## Calculate h,P
        self.h = math.sqrt(MEW * self.a * (1 - math.pow(self.e, 2)))
        self.P = math.pow(self.h, 2) / MEW

        ## Initialize sgp4 object
        self.sgpObject = twoline2rv(line1, line2, wgs72)

        self.a_s.append(self.a)
        self.i_s.append(self.i)

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
            raise PropagationError(self.sgpObject.error_message, self.sgpObject.error, time)

        return np.array(r_teme)[np.newaxis], np.array(v_teme)[np.newaxis]

    def get_eci_state(self, time, matlabLock):
        """
        At the specified time, return the state vector in ECI frame

        :param datetime.datetime time: time of state vector
        :param threading.RLock matlabLock: lock for matlab access
        :return: numpy.array r_eci, numpy.array v_eci
        """
        r_teme, v_teme = self.get_teme_state(time)

        r_eci, v_eci = teme2eci(r_teme.T, v_teme.T, time, delta_at, matlabLock)

        return r_eci.T, v_eci.T

    def generate_trajectory(self, startTime, endTime, steps, matlabLock, laserObject = None):
        """
        Generate a trajectory between startTime and endTime in interval time steps

        :param datetime.datetime startTime: start of trajectory (inclusive)
        :param datetime.datetime endTime: end of trajectory (inclusive)
        :param int steps: number of time steps between trajectory points
        :param threading.RLock matlabLock: lock to ensure one thread calls MATLAB at a time
        :param Laser laserObject: laser object in orbit OPTIONAL
        :return: define pandas.DataFrame self.trajectory
        """

        deltaTime = int((endTime - startTime).total_seconds())
        self.timeStep = deltaTime/steps
        step = round(deltaTime / steps)
        times = list(range(0, deltaTime, step))
        times.append(deltaTime)

        self.trajectoryTimes = []
        self.trajectoryPos = np.zeros([len(times), 3])
        self.trajectoryVeloc = np.zeros([len(times), 3])

        for indx, deltaSeconds in enumerate(times):
            time = startTime + dt.timedelta(seconds=deltaSeconds)

            r, v = self.get_teme_state(time)

            if laserObject is not None and laserObject.enable:
                laserPos, laserVel = laserObject.parse_trajectory(time=time)
                relPos = (r - np.array(laserPos))[0]
                relDist = np.linalg.norm(relPos)
                fireDuration = 60 #sec

                if relDist < laserObject.range and laserObject.is_ready(time, fireDuration):
                    print('Firing laser on {} ({}) at {}'.format(self.satName, self.satNum, time.strftime('%Y-%m-%d %H:%M:%S')))
                    deltaV = laserObject.fire(time, time+dt.timedelta(seconds=fireDuration), self)
                    unitPos = relPos / relDist
                    v = v + deltaV*unitPos

                    r_ECI, v_ECI = teme2eci(r.T, v.T, time, delta_at, matlabLock)

                    self.tle = self.update_tle(r_ECI.T, v_ECI.T, matlabLock, time)
                    self.parse_tle()

            self.trajectoryTimes.append(time)
            self.trajectoryPos[indx,:] = r
            self.trajectoryVeloc[indx,:] = v

        combined_rv = np.hstack((self.trajectoryPos, self.trajectoryVeloc))
        self.trajectory = pd.DataFrame(data=combined_rv, columns=['Posx','Posy','Posz','Velx','Vely','Velz'])
        self.trajectory['Times'] = self.trajectoryTimes

    def update_tle(self, r, v, matlabLock, time):
        """
        Return an updated TLE from STK using the input state vectors in ECI

        :param numpy.array r: radius in ECI [m]
        :param numpy.array v: velocity in ECI [m]
        :param threading.RLock matlabLock: lock for matlab access
        :param datetime.datetime time: current time
        :return: string tle: new TLE
        """

        ## Convert to km
        r = r/1000
        v = v/1000
        arry = [time.year, time.month, time.day, 17, 0, 0]

        with matlabLock:
            line1, line2 = eng.tle_from_stk(matlab.double(r.tolist()), matlab.double(v.tolist()),
                                           matlab.int64([self.satNum]), matlab.double(arry), nargout=2)
        tle = '0 {}\n{}\n{}'.format(self.satName, line1, line2)
        return tle

    def parse_trajectory(self, indx = None, time = None):
        """
        Given index or time in the trajectory DataFrame, return r, v in nice form.

        :param int indx: indx in self._trajectory OPTIONAL
        :param datetime.datetime tim: time in self._trajectory OPTIONAL
        :return: float list r, float list v in format [x,y,z]
        """

        if indx is None and time is None:
            raise ValueError("Either indx or time argument must be given")

        if indx is not None and time is not None:
            raise ValueError("Either indx argument or time argument must be none, both cannot be given")

        r = []
        v = []
        if time is not None:
            ##TODO: Generalize for position closest to given time (if time is not in trajectory)
            entry = self.trajectory.loc[self.trajectory['Times'] == time]

            if len(entry) == 0:
                return r, v

            row = entry.iloc[0]
            r = [row['Posx'], row['Posy'],
                 row['Posz']]
            v = [row['Velx'], row['Vely'],
                 row['Velz']]

        if indx is not None:
            r = [self.trajectory.iloc[indx]['Posx'], self.trajectory.iloc[indx]['Posy'],
                 self.trajectory.iloc[indx]['Posz']]
            v = [self.trajectory.iloc[indx]['Velx'], self.trajectory.iloc[indx]['Vely'],
                 self.trajectory.iloc[indx]['Velz']]

        return r, v

class Laser(Object):
    """
    Object class specifically for a Laser
    """
    chargingRate = 380.75 #W
    laserPower = 32.3e3 #kW
    range = 100e3
    deltaV_persec = 0.4  # m/s
    def __init__(self, obj, enable):
        """
        Init class

        :param Object obj: input obj that defines parameters
        :param bool enable: enable laser in simulation
        """

        super().__init__(obj.tle)
        self.enable = enable

        self._fireTimes = pd.DataFrame([],columns=['Start','End','Duration'])
        self._fireObjs = []

    def fire(self, startTime, endTime, deb):
        """
        Return the delta-V generated by a laser firing on the object

        :param datetime.datetime startTime: beginning time of laser fire
        :param datetime.datetime endTime: ending time of laser fire
        :param float range: distance between laser and debris
        :param Object deb: Object of debris fired on
        :return: float deltaV: deltaV generated by laser on object [m/s]
        """

        tempS = pd.Series([startTime, endTime, (endTime-startTime).total_seconds()],index=['Start','End','Duration'])
        self._fireTimes = append(self._fireTimes, tempS)
        self._fireTimes = self._fireTimes.sort_values(by='Start')

        if deb not in self._fireObjs:
            self._fireObjs.append(deb)

        #TODO: Account for angle of laser -> don't fire if bad angle
        deltaV = self.deltaV_persec * (endTime - startTime).total_seconds()

        return deltaV

    def is_ready(self, time, duration):
        """
        Return whether the laser is charged and able to fire at time

        :param datetime.datetime time: time of potential firing
        :param float duration: length of proposed laser firing [sec]
        :return: bool ready
        """
        if not self.enable:
            return False

        if len(self._fireTimes) != 0:
            relTimes = np.array([(self._fireTimes.iloc[i]['Start'] - time).total_seconds()
                        for i in range(0,len(self._fireTimes))])
            indices = np.where(relTimes > 0)

            if len(indices[0]) == 0:
                indices = np.where(relTimes < 0)

                if len(indices[0]) == 0:
                    return False

                changeOver = indices[0][-1]
                beforeRow = self._fireTimes.iloc[changeOver]
                beforeTime = beforeRow['End']
                delta = (beforeTime - time).total_seconds()
                nrgUsed = self.laserPower * beforeRow['Duration']
            else:
                changeOver = indices[0][0]

                afterTime = self._fireTimes.iloc[changeOver]['Start']
                if changeOver >= 1:
                    beforeRow = self._fireTimes.iloc[changeOver-1]
                    beforeTime = beforeRow['End']

                    delta1 = (afterTime - time).total_seconds()
                    nrgUsed1 = self.laserPower*duration

                    delta2 = (time - beforeTime).total_seconds()
                    nrgUsed2 = self.laserPower*beforeRow['Duration']

                    return self.ready(nrgUsed1, delta1) and self.ready(nrgUsed2, delta2)

                else:
                    delta = (afterTime - time).total_seconds()
                    nrgUsed = self.laserPower*duration

            return self.ready(nrgUsed, delta)

        return True

    def ready(self, nrgUsed, delta):
        """
        Given energy used and the elapsed time, calculate if the energy has been re-charged

        :param float nrgUsed: energy used in the last pulse
        :param float delta: elapsed time in seconds
        :return: bool ready: True if laser is ready to fire again
        """
        ## TODO: Better account for time in eclipse
        chargeTime = nrgUsed / self.chargingRate * 2 #multiply by 2 to account for time in eclipse, on average
        if chargeTime > abs(delta):
            return False

        return True

    def get_fire_times(self):
        """
        Return fireTimes dataframe

        :return: pandas.DataFrame tracking fire times
        """
        return self._fireTimes

    def get_fire_objs(self):
        """
        Return fireObjs list
        :return: list int fireObjs
        """
        return self._fireObjs

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

def get_jd(time):
    """
    Return julian date at given time

    :param datetime.datetime time: time at which to calculate JD
    :return: float JD
    """

    ##TODO: Account for leapseconds
    yr = time.year
    month = time.month
    day = time.day
    hour = time.hour
    min = time.minute
    sec = time.second

    JD = 367*yr - int((7*(yr+int((month+9)/12))/4)) + int(275*month/9) + day + 1721013.5 + ((sec/60+min)/60+hour)/24
    return JD

def teme2eci(r_teme, v_teme, time, delta_at, matlabLock):
    """
    Convert input r and v vectors to ECI frame (j2000)

    :param numpy.array r_teme: radius in teme, col vector format
    :param numpy.array v_teme: velocity in teme, col vector format
    :param datetime.datetime time: time of conversion
    :param datetime.timedelta delta_at: difference from TAI to UT1
    :param threading.RLock matlabLock: lock to access matlab from one thread
    :return: numpy.array r, numpy.array v in ECI frame
    """
    ## Get epoch time
    tai = time + delta_at
    tt = tai + dt.timedelta(seconds=32.184)

    ttt = (get_jd(tt) - 2451545) / 36525
    with matlabLock:
        r_eci, v_eci, aeci = eng.teme2eci(matlab.double(r_teme.tolist()), matlab.double(v_teme.tolist()),
                                          matlab.double([[0], [0], [0]]), ttt, matlab.double([0]),
                                          matlab.double([0]), nargout=3)

    return np.asarray(r_eci), np.asarray(v_eci)

def row2col(vec):
    """
    Convert 3 element vector from row into column

    :param numpy.array vec: row vector to rotate
    :return: list col
    """
    col = [[vec[0]], [vec[1]], [vec[2]]]

    return col

def rnd(num, delim):
    """
    Rounding function to round number to nearest even value divisible by delim

    :param float num: number to round
    :param float delim: difference between even values
    :return: int rounded
    """

    lower_even = int(num // delim) * delim
    upper_even = lower_even + delim

    if abs(num-lower_even) < abs(num-upper_even):
        return int(lower_even)

    return int(upper_even)

def get_bounds(arry, x_arry, y_arry):
    """
    Given arry, return the bounds for x and y to set the plot (i.e. find the lowest and highest non-zero values of x
    and y in the data).

    :param numpy.array arry: data array
    :param numpy.array x_arry: range of x values
    :param numpy.array y_arry: range of y values
    :return: tuple x_bounds, tuple y_bounds in format (min_val, max_val)
    """
    bounds = arry.nonzero()
    bounds[0].sort()
    bounds[1].sort()
    if len(bounds[0]) == 0:
        min_y = 0
        max_y = -1
    else:
        min_y = bounds[0][0]
        max_y = bounds[0][-1]

    if len(bounds[1]) == 0:
        min_x = 0
        max_x = -1
    else:
        min_x = bounds[1][0]
        max_x = bounds[1][-1]

    if min_x == max_x:
        min_x = 0
        max_x = -1

    if min_y == max_y:
        min_y = 0
        max_y = -1

    return (x_arry[min_x], x_arry[max_x]), (y_arry[min_y], y_arry[max_y])

def append(exist, new):
    """
    Append dataframe to end of existing dataframe

    :param pd.DataFrame exist: exisiting 'large' df
    :param pd.DataFrame or pd.Series new: new 'small' df or series to append to end of existing
    :return: pd.DataFrame appended, combination of the two
    """

    if type(new) == pd.DataFrame:
        appended = pd.concat([exist, new], axis=0, ignore_index=True)
    else:
        appended = pd.concat([exist.T, new], axis=1, ignore_index=True)
        return appended.T

    return appended

def split_array(tuple):
    """
    Split array into sub-arrays with continuous values

    :param tuple array: input tuple with first element an np.array to split
    :return: np.array split
    """
    split = np.split(tuple[0], np.where(np.diff(tuple[0]) != 1)[0] + 1)

    for indx, array in enumerate(split):
        if array.ndim == 2 and len(array) > 0:
            split[indx] = split[indx][0]

    return split