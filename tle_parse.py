import math
import datetime as dt
import argparse
import pickle
import os

import numpy as np
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool, Legend, LogColorMapper, ColorBar, LogTicker
from bokeh.layouts import gridplot

class Object:
    """
    Generic class for orbital object
    """
    def __init__(self, tle_str):
        """
        Initializer for object class
        :param tle_str: tle in string format, used to define object params
        """

        lines = tle_str.split('\n')
        line1 = lines[0]
        line2 = lines[1]

        self.sat_num = int(line1[2:7])
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
        :return: Defne self.trajectory
        """

        self._T = 2*math.pi/self.n
        endTime = int(self._T//self._intervalTime*self._intervalTime) #round T to nearest interval

        times = range(0, endTime+self._intervalTime, self._intervalTime)
        self._trajectory = np.zeros([len(times),3])
        for i in range(0,len(times)):
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

    def is_in_sphere(self,laser_loc, time, sphere_size):
        """
        Given the location of the laser satellite, the time, and sphere size, calculates if obj is in sphere
        :param list laser_loc: location of laser in ECI [x,y,z] [m]
        :param datetime time: time in UTC
        :param float sphere_size: radius of targeting sphere [m]
        :return: bool in: True if obj is in the sphere
        """

        ## Find object location at time
        allTimeSincePerigee = abs((time - self._pt).total_seconds()) #absolute time since last perigee passage
        orbitTimeSincePerigee = allTimeSincePerigee % self._T #time within 1 orbit since perigee passage
        closestTimeIndex = int(orbitTimeSincePerigee//self._intervalTime)
        closestPoint = self._trajectory[closestTimeIndex,:]

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

def parse_catalog():
    """
    Open the txt file with all TLEs in the sat cat, create an Object for each tle
    :return: objects list, each index is an Object in the sat cat
    """
    with open('satcat.txt','r') as f:
        lines = f.readlines()
        objects = [None]*int(len(lines)/2)
        for i in range(0, len(lines), 2):
            tle = lines[i:i+2]
            str = '{}{}'.format(tle[0], tle[1])
            objects[int(i/2)] = Object(str)

    return objects

MEW = 3986004.418e8 #m3/s2
Re = 6378.37e3 # m

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Input arguments")
    required = parser.add_argument_group('optional arguments')
    required.add_argument("--parse_data", "-p", help="Re-parse sat cat file, rather than read stored data",
                          required=False, action="store_true", default=False)
    #optional = parser.add_argument_group('optional arguments')
    args = vars(parser.parse_args())

    parse = args['parse_data']
    savefile = 'satcat.pickle'
    if parse or savefile not in os.listdir(os.getcwd()):
        print('Parsing tle file...')
        objects = parse_catalog()
        with open(savefile,'wb') as f:
            pickle.dump(objects,f)
    else:
        print('Loading pickle file')
        with open(savefile,'rb') as f:
            objects = pickle.load(f)

    ## Check that Plots folder is made
    if not 'Plots' in os.listdir(os.getcwd()):
        os.makedirs(os.getcwd()+'/Plots')

    lower_alt = 200
    upper_alt = 40e3
    alt_delim = 50
    lower_i = 0
    upper_i = 180
    i_delim = 5

    np_alt = np.linspace(lower_alt, upper_alt, int((upper_alt+alt_delim-lower_alt)/alt_delim))
    altitudes = dict.fromkeys(list(np_alt),0)
    np_inc = np.linspace(lower_i, upper_i, int((upper_i+i_delim-lower_i)/i_delim))
    inclinations = dict.fromkeys(list(np_inc),0) #Range from 0 to 180 deg inclin. in 5 deg steps

    z = np.zeros([len(np_inc),len(np_alt)])

    for obj in objects:
        inc = math.degrees(obj.i)
        alt = (obj.a - Re)/1000 #km

        nearest_alt = int(alt//alt_delim)*alt_delim
        nearest_inc = int(inc//i_delim)*i_delim

        if nearest_alt >= lower_alt and nearest_alt <= upper_alt:
            altitudes[nearest_alt]+= 1

        if nearest_inc >= lower_i and nearest_inc <= upper_i:
            inclinations[nearest_inc]+=1

        ## We need both within range to generate the image array
        if nearest_alt < lower_alt or nearest_alt > upper_alt \
            or nearest_inc < lower_i or nearest_inc > upper_i:
            continue

        x_coord = np.where(np_alt == nearest_alt)[0][0]
        y_coord = np.where(np_inc == nearest_inc)[0][0]

        z[y_coord][x_coord]+= 1


    cmap = LogColorMapper(palette='Viridis256',low=1,high=z.max())
    output_file('Plots/tle_distro.html')
    p = figure(title='Distribution of TLEs', x_axis_label='Altitude (km)', y_axis_label='Debris Amount')
    p.vbar(x=list(altitudes.keys()), top=list(altitudes.values()), bottom=0, width=alt_delim)

    p2 = figure(title='Distribution of TLEs', x_axis_label='Inclination (deg)', y_axis_label='Debris Amount')
    p2.vbar(x=list(inclinations.keys()), top=list(inclinations.values()), bottom=0, width=i_delim)

    p3 = figure(title='Distribution of TLEs',x_axis_label='Altitude (km)',y_axis_label='Inclination (deg)',
                x_range=(lower_alt,upper_alt),y_range=(lower_i,upper_i),
                tooltips=[("x","$x"),("y","$y"),("value","@image")])
    p3.image(image=[z],x=lower_alt,y=lower_i,dw=upper_alt-lower_alt,dh=upper_i-lower_i,color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap,location=(0,0),ticker=LogTicker())
    p3.add_layout(color_bar,'right')
    grid = gridplot([[p, p2],[p3]])
    show(grid)