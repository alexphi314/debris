## Alex Philpott
## Run an orbit propagator for all objects and the laser.
## Generates almost all the plots used for orbit determination.

import argparse
import pickle
import os
import datetime as dt
import threading
import queue
import sys
import math

import pandas as pd
import numpy as np
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.palettes import Blues9
from bokeh.layouts import gridplot

from astroUtils import parse_catalog, Re, PropagationError, Object, Laser, rnd, get_bounds, append, split_array

frmat = '%Y-%m-%d %H:%M:%S'

class Simulator:
    """
    Simulator manager class
    """
    def __init__(self, startTime, endTime, steps):
        """
        :param datetime.datetime startTime: simulation start time
        :param datetime.datetime endTime: simulation end time
        :param int steps: number of trajectory steps
        """

        ## Define trajectory variables
        self.trajQueue = queue.Queue()
        self.startTime = startTime
        self.endTime = endTime
        self.steps = steps
        self.badTrajs = []
        self.dfCols = ['Time', 'Object', 'Distance [km]', 'Number', 'Rel V', 'Deb I']
        self.decayed = pd.DataFrame([],columns=['Time','Number'])

        self.passQueue = queue.Queue()
        self.aLock = threading.Lock()
        self.durationALock = threading.Lock()
        self.matlabLock = threading.RLock()
        self.decayedALock = threading.Lock()

    def setLaser(self, laser):
        """
        Update the laser object used for distance calcs

        :param Laser obj: sat cat obj used as laser
        :return:
        """
        assert self.passQueue.empty()
        assert self.trajQueue.empty()


        self.laserObject = laser
        self.laserObject.generate_trajectory(self.startTime, self.endTime, self.steps, self.matlabLock)

        ## Define laser pass variables
        self._laserPasses = pd.DataFrame([], columns=self.dfCols)
        self._passDurations = pd.DataFrame([], columns=['Start Time','End Time','Duration'])

    def gen_trajectories(self):
        """
        For all objects in queue, generate trajectory

        :return:
        """
        while True:
            try:
                obj = self.trajQueue.get()
                obj.generate_trajectory(self.startTime, self.endTime, self.steps, self.matlabLock, self.laserObject)
                #self.message('Trajectory generated for {}'.format(obj.satName))
            except PropagationError as e:
                self.message('Got Propagation Error {}: {}'.format(e.num, e.msg))
                self.badTrajs.append(obj)

                if e.num == 6:
                    if len(self.decayed) == 0:
                        num = 1
                    else:
                        num = self.decayed.iloc[-1]['Number']+1
                    time = (e.time-self.startTime).total_seconds()
                    temp = pd.Series([time/3600/24,num],index=['Time','Number'])
                    with self.decayedALock:
                        self.decayed = append(self.decayed, temp)
                pass
            self.trajQueue.task_done()

    def compute_passes(self):
        """
        For all objects in queue, calculate distance to laser at each time step

        :return:
        """
        while True:
            try:
                obj = self.passQueue.get()
                assert self.laserObject is not None

                ## Do not include passes of non-debris
                if not obj.deb:
                    self.passQueue.task_done()
                    continue

                relPosition = obj.trajectoryPos - self.laserObject.trajectoryPos
                relVelocity = obj.trajectoryVeloc - self.laserObject.trajectoryVeloc

                relDist = np.linalg.norm(relPosition, axis=1)
                relVel = np.linalg.norm(relVelocity, axis=1)
                indices = np.where(relDist < self.laserObject.range)
                smallDist = relDist[indices]/1000 #convert to km
                smallVel = relVel[indices]/1000 #conver to km/s

                if len(smallDist) > 0:
                    tempDf = pd.DataFrame(columns=self.dfCols)
                    tempDf['Distance [km]'] = smallDist
                    tempDf['Time'] = [obj.trajectoryTimes[i] for i in indices[0]]
                    tempDf['Object'] = obj.satName
                    tempDf['Number'] = obj.satNum
                    tempDf['Rel V'] = smallVel
                    tempDf['Deb I'] = math.degrees(obj.i)

                    if len(smallDist) > 1:
                        split = split_array(indices)

                        for arry in split:
                            if len(arry) == 1:
                                continue

                            accessStart = obj.trajectoryTimes[arry[0]]
                            accessEnd = obj.trajectoryTimes[arry[-1]]

                            tempSeries = pd.Series([accessStart, accessEnd, (accessEnd - accessStart).total_seconds()],
                                                   index=['Start Time', 'End Time', 'Duration'])

                            with self.durationALock:
                                self._passDurations = append(self._passDurations, tempSeries)

                    assert len([obj.trajectoryTimes[i] for i in indices[0]]) == len(smallDist)

                    with self.aLock:
                        self._laserPasses = append(self._laserPasses, tempDf)

            except Exception as e:
                e_type, e_obj, e_tb = sys.exc_info()
                self.message('Got {} on line {}: {}'.format(e_type, e_tb.tb_lineno, e))
                pass
            self.passQueue.task_done()

    def get_passes(self):
        """
        Return laserPasses dataframe

        :return: pd.DataFrame laserPasses
        """
        return self._laserPasses

    def get_durations(self):
        """
        Return passDurations dataframe

        :return: pd.DataFrame passDurations
        """
        return self._passDurations

    def message(self, msg):
        """
        Print message with thread name

        :param msg: message to print
        :return:
        """
        print('{}: {}'.format(threading.current_thread().name, msg))

if __name__ == "__main__":

    #####################################
    ### Read Arguments and Load Files ###
    #####################################

    parser = argparse.ArgumentParser(description="Input arguments")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--parse_data", "-p", help="Re-parse sat cat file, rather than read stored data",
                          required=False, action="store_true", default=False)
    optional.add_argument("--sweep_objs", "-s", help="Sweep through satellites for potential lasers",
                          required=False, action="store_true", default=False)
    args = vars(parser.parse_args())

    print('Running with call: {}'.format(sys.argv[0:]))

    parse = args['parse_data']
    sweep = args['sweep_objs']
    dataDir = 'Data'
    savefile1 = 'satcat3_LEO_deb.pickle'
    savefile2 = 'satcat3_LEO_obj.pickle'
    if parse or savefile1 not in os.listdir(os.path.join(os.getcwd(),dataDir))\
            or savefile2 not in os.listdir(os.path.join(os.getcwd(),dataDir)):
        print('Parsing tle file...')
        objects = parse_catalog(3)

        ## Throw out non-LEO objects
        deb = []
        sats = []
        for i in range(0, len(objects)):
            obj = objects[i]
            h = (obj.a - Re) / 1e3  # km

            if h < 2000:
                if obj.deb:
                    deb.append(obj)
                else:
                    sats.append(obj)

        with open(os.path.join(os.getcwd(),dataDir,savefile1), 'wb') as f:
            pickle.dump(deb, f)
        with open(os.path.join(os.getcwd(),dataDir,savefile2), 'wb') as f:
            pickle.dump(sats, f)
    else:
        print('Loading pickle files')
        with open(os.path.join(os.getcwd(),dataDir,savefile1), 'rb') as f:
            deb = pickle.load(f)
        with open(os.path.join(os.getcwd(), dataDir, savefile2), 'rb') as f:
            sats = pickle.load(f)

    ## Check that folders are made
    if not 'Plots' in os.listdir(os.getcwd()):
        os.makedirs(os.getcwd() + '/Plots')

    if not 'Data' in os.listdir(os.getcwd()):
        os.makedirs(os.getcwd() + '/Data')

    ######################
    ### Run Simulation ###
    ######################

    print('Starting sim...')
    print('Running with {} pieces of debris and {} sats'.format(len(deb), len(sats)))
    objects = deb + sats
    numDays = 7
    startTime = dt.datetime(2019,3,27,17,00,00)
    endTime = startTime + dt.timedelta(days=numDays)
    steps = numDays*6*24
    enableLaser = False

    laser = None
    laserIndx = None

    # Target Laser Object: METOP-C
    laser_objs = []
    for indx,obj in enumerate(objects):
        if obj.satNum == 43689:
            laser = Laser(obj, enableLaser)
            laserIndx = indx
            laser_objs.append(laser)
            break
    del objects[laserIndx]

    simulator = Simulator(startTime, endTime, steps)
    simulator.setLaser(laser)

    for i in range(1,5):
        worker = threading.Thread(target=simulator.gen_trajectories,
                                  name='trajectory-worker-{}'.format(i))
        worker.setDaemon(True)
        worker.start()
    for i in range(1,5):
        worker = threading.Thread(target=simulator.compute_passes,
                                  name='pass-worker-{}'.format(i))
        worker.setDaemon(True)
        worker.start()

    for obj in objects:
        simulator.trajQueue.put(obj)

    simulator.trajQueue.join()
    print('Finished generating trajectories, with {} invalid'.format(
        len(simulator.badTrajs))
    )

    ## Generate range arrays
    lower_alt = 695
    upper_alt = 1005
    alt_delim = 1
    np_alt = np.linspace(lower_alt, upper_alt, int((upper_alt - lower_alt) / alt_delim) + 1)

    lower_i = -1
    upper_i = 181 #boundaries for actual debris objects
    i_delim = 1
    np_inc = np.linspace(lower_i, upper_i, int((upper_i - lower_i) / i_delim) + 1)

    lower_v = 0
    upper_v = 25
    v_delim = 0.25
    np_vel = np.linspace(lower_v, upper_v, int((upper_v - lower_v) / v_delim) + 1)

    slower_v = 0
    supper_v = 50/1000
    sv_delim = 0.001
    np_svel = np.linspace(slower_v, supper_v, int((supper_v - slower_v) / sv_delim) + 1)

    lower_d = 0
    upper_d = 100
    d_delim = 1
    np_dist = np.linspace(lower_d, upper_d, int((upper_d - lower_d) / d_delim) + 1)

    uniq = np.zeros([len(np_inc), len(np_alt)])
    tot = uniq
    vel = np.zeros([len(np_vel), len(np_dist)])
    inc_arry = np.zeros([len(np_inc), len(np_inc)])
    svel = np.zeros([len(np_svel), len(np_dist)])

    if sweep:
        i_sat_bounds = [97, 100]  # boundary inclinations for laser objects
        alt_sat_bounds = [810, 830]
        np_alt_sat = np.linspace(alt_sat_bounds[0], alt_sat_bounds[1],
                                 int((alt_sat_bounds[1] - alt_sat_bounds[0]) / 1) + 1)
        np_inc_sat = np.linspace(i_sat_bounds[0], i_sat_bounds[1], int((upper_i - lower_i) / i_delim) + 1)

        appended = []
        for alt in np_alt_sat:
            for inc in np_inc_sat:
                for obj in objects:
                    if 'TBA' in obj.satName or obj.deb:
                        continue

                    h = (obj.a - Re) / 1000
                    i = math.degrees(obj.i)

                    if abs(h-alt) <= alt_delim/2 and abs(i-inc) <= i_delim/2 and obj not in appended:
                        # print('{} ({}): alt {} km inc {} deg ecc {}'.format(
                        #     obj.satName, obj.satNum, round(h, 1), round(i, 1), round(obj.e, 2)
                        # ))
                        laser_objs.append(Laser(obj.tle, enableLaser))
                        appended.append(obj)
                        break

    else:
        target_a = 820
        target_i = 99
        a_tol = 1
        i_tol = 0.5
        # for obj in objects:
        #     if 'TBA' in obj.satName or obj.deb:
        #         continue
        #
        #     h = (obj.a - Re) / 1000
        #     i = math.degrees(obj.i)
        #
        #     if abs(h - target_a) <= a_tol and abs(i - target_i) <= i_tol and obj not in laser_objs:
        #         # print('{} ({}): alt {} km inc {} deg ecc {}'.format(
        #         #     obj.satName, obj.satNum, round(h, 1), round(i, 1), round(obj.e, 2)
        #         # ))
        #         laser_objs.append(obj)

        # Target Laser Object: METOP-C
        # for obj in objects:
        #     if obj.satNum == 43689:
        #         laser_objs.append(Laser(obj.tle, enableLaser))
        #         break

    Blues9.reverse()

    print('Running from {} to {}, with {} steps'.format(
        startTime.strftime(frmat), endTime.strftime(frmat), steps
    ))

    for obj in objects:
        if obj not in simulator.badTrajs:
            simulator.passQueue.put(obj)

    simulator.passQueue.join()

    ###############################
    ### Output Results and Plot ###
    ###############################

    passes = simulator.get_passes()

    durations = simulator.get_durations()
    lower_dur = 0
    upper_dur = durations['Duration'].max()
    dur_delim = 60
    np_dur = np.linspace(lower_dur, upper_dur, int((upper_dur - lower_dur) / dur_delim) + 1)
    durs = np.zeros([len(np_dur)])

    print('')
    print('Laser Object: {} ({})'.format(simulator.laserObject.satName, simulator.laserObject.satNum))
    uniquePasses = set([passes.iloc[i]['Number'] for i in range(0,len(passes))])
    print('Got {} passes, {} unique objects'.format(len(passes), len(uniquePasses)))

    alt = (simulator.laserObject.a - Re)/1000
    nearest_alt = rnd(alt, alt_delim)
    nearest_inc = rnd(math.degrees(simulator.laserObject.i), i_delim)

    x_coord = np.where(np_alt == nearest_alt)[0][0]
    y_coord = np.where(np_inc == nearest_inc)[0][0]
    uniq[y_coord][x_coord] += len(uniquePasses)
    tot[y_coord][x_coord] += len(passes)

    for indx, row in passes.iterrows():
        nearest_dist = rnd(row['Distance [km]'], d_delim)
        nearest_deb_i = rnd(row['Deb I'], i_delim)
        nearest_vel = rnd(row['Rel V'], v_delim)

        dist_coord = np.where(np_dist == nearest_dist)[0][0]
        deb_inc_coord = np.where(np_inc == nearest_deb_i)[0][0]
        vel_coord = np.where(np_vel == nearest_vel)[0][0]

        if row['Rel V'] < supper_v:
            nearest_svel = rnd(row['Rel V'], sv_delim)
            svel_coord = np.where(np_svel == nearest_svel)[0][0]
            svel[svel_coord][dist_coord] += 1

        vel[vel_coord][dist_coord]+=1
        inc_arry[deb_inc_coord][y_coord]+=1

    for indx, row in durations.iterrows():
        nearest_dur = rnd(row['Duration'], dur_delim)
        coord = np.where(np_dur == nearest_dur)[0][0]
        durs[coord] += 1

    output_file('Plots/passDurations_{}_{}.html'.format(simulator.laserObject.satNum,steps))
    p6 =figure(title='Pass Duration Distribution', x_axis_label='Duration (sec)', y_axis_label='Number of Passes',
                x_range=(lower_dur, upper_dur), y_range=(0, durs.max()),
                tooltips=[("Duration", "$x"), ("Number", "$y")])
    p6.vbar(x=np_dur, top=durs, bottom=0, width=dur_delim)

    delta_as = []
    delta_is = []

    for obj in simulator.laserObject.get_fire_objs():
        for i in range(1,len(obj.a_s)):
            deltaA = obj.a_s[i] - obj.a_s[i-1]
            delta_as.append(deltaA)

        for i in range(1,len(obj.i_s)):
            deltaI = obj.i_s[i] - obj.i_s[i-1]
            delta_is.append(math.degrees(deltaI))

    p9 = figure(title='Orbit Changes', x_axis_label='Change in Semi-Major Axis [m]',
                y_axis_label='Change in Inclination [deg]')
    p9.scatter(x=delta_as,y=delta_is)

    grid = gridplot([[p6,p9]])
    show(grid)

    uniq_x, uniq_y = get_bounds(uniq, np_alt, np_inc)
    tot_x, tot_y = get_bounds(tot, np_alt, np_inc)
    vel_x, vel_y = get_bounds(vel, np_dist, np_vel)
    inc_x, inc_y = get_bounds(inc_arry, np_inc, np_inc)
    svel_x, svel_y = get_bounds(svel, np_dist, np_svel)

    print('Uniq: x ({} {}) y ({} {})'.format(uniq_x[0], uniq_x[1], uniq_y[0], uniq_y[1]))
    print('Tot: x ({} {}) y ({} {})'.format(tot_x[0], tot_x[1], tot_y[0], tot_y[1]))
    print('Vel: x ({} {}) y ({} {})'.format(vel_x[0], vel_x[1], vel_y[0], vel_y[1]))
    print('Inc: x ({} {}) y ({} {})'.format(inc_x[0], inc_x[1], inc_y[0], inc_y[1]))

    output_file('Plots/laserPasses_{}.html'.format(steps))
    p2 = figure(title='Debris Passes', x_axis_label='Altitude (km)', y_axis_label='Inclination (deg)',
                x_range=(uniq_x[0], uniq_x[1]), y_range=(uniq_y[0], uniq_y[1]),
                tooltips=[("Alt", "$x"), ("Inc", "$y"), ("Passes", "@image")])
    cmap = LinearColorMapper(palette=Blues9, low=0, high=uniq.max())
    p2.image(image=[uniq], x=lower_alt, y=lower_i, dw=upper_alt - lower_alt, dh=upper_i - lower_i, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0))
    p2.add_layout(color_bar, 'right')

    p3 = figure(title='Total Passes', x_axis_label='Altitude (km)', y_axis_label='Inclination (deg)',
                x_range=(tot_x[0], tot_x[1]), y_range=(tot_y[0], tot_y[1]),
                tooltips=[("Alt", "$x"), ("Inc", "$y"), ("Passes", "@image")])
    cmap = LinearColorMapper(palette=Blues9, low=0, high=tot.max())
    p3.image(image=[tot], x=lower_alt, y=lower_i, dw=upper_alt - lower_alt, dh=upper_i - lower_i, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0))
    p3.add_layout(color_bar, 'right')

    p4 = figure(title='Rel V vs. Distance', x_axis_label='Distance (km)', y_axis_label='Relative Velocity (km/s)',
                x_range=(vel_x[0], vel_x[1]), y_range=(vel_y[0], vel_y[1]),
                tooltips=[("Dist", "$x"), ("Vel", "$y"), ("Number", "@image")])
    cmap = LinearColorMapper(palette=Blues9, low=0, high=vel.max())
    p4.image(image=[vel], x=lower_d, y=lower_v, dw=upper_d-lower_d, dh=upper_v-lower_v, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0))
    p4.add_layout(color_bar, 'right')

    p5 = figure(title='Obj Inc vs. Laser Inc', x_axis_label='Laser Inclination (deg)', y_axis_label='Debris Inclination (deg)',
                x_range=(inc_x[0], inc_x[1]), y_range=(inc_y[0], inc_y[1]),
                tooltips=[("Laser Inc", "$x"), ("Deb Inc", "$y"), ("Number", "@image")])
    cmap = LinearColorMapper(palette=Blues9, low=0, high=inc_arry.max())
    p5.image(image=[inc_arry], x=lower_i, y=lower_i, dw=upper_i-lower_i, dh=upper_i-lower_i, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0))
    p5.add_layout(color_bar, 'right')

    source = ColumnDataSource(simulator.decayed)
    p7 = figure(title='Debris Decays', x_axis_label='Elapsed Time [days]',
                y_axis_label='Number of Debris Decays',tooltips=[('Time','$x'),('Number','$y')])
    p7.line(x=list(simulator.decayed['Time']),y=list(simulator.decayed['Number']))

    p8 = figure(title='Zoomed Relative V vs. Distance', x_axis_label='Distance (km)', y_axis_label='Relative Velocity (km/s)',
                x_range=(svel_x[0], svel_x[1]), y_range=(svel_y[0], svel_y[1]),
                tooltips=[("Dist", "$x"), ("Vel", "$y"), ("Number", "@image")])
    cmap = LinearColorMapper(palette=Blues9, low=0, high=svel.max())
    p8.image(image=[svel], x=lower_d, y=slower_v, dw=upper_d - lower_d, dh=supper_v - slower_v, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0))
    p8.add_layout(color_bar, 'right')

    grid = gridplot([[p2,p3],[p4,p5],[p7,p8]])
    show(grid)