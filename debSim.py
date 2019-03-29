## Run collision simulator

import argparse
import pickle
import os
import datetime as dt
import threading
import queue
import sys

import pandas as pd
import numpy as np
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource

from astroUtils import parse_catalog, Re, PropagationError, Object

frmat = '%Y-%m-%d %H:%M:%S'

def append(exist,new):
    """
    Append dataframe to end of existing dataframe
    :param pd.DataFrame exist: exisiting 'large' df
    :param pd.DataFrame new: new 'small' df to append to end of existing
    :return: pd.DataFrame appended, combination of the two
    """
    if len(exist) == 0:
        return new

    appended = pd.concat([exist,new],axis=0,ignore_index=True)
    return appended

class Simulator:
    """
    Simulator manager class
    """
    tempLaserTLE = "LASER\n{}\n{}".format(
        '1   900U 64063C   19085.91254635 +.00000212 +00000-0 +21859-3 0  9999',
        '2   900 090.1517 021.2926 0026700 330.0700 101.6591 13.73206334708891'
    )
    laserRange = 100e3 #m
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
        self.dfCols = ['Time', 'Object', 'Distance [km]', 'Number']

        self.passQueue = queue.Queue()
        self.aLock = threading.Lock()

    def setLaser(self, obj):
        """
        Update the laser object used for distance calcs
        :param Object obj: sat cat obj used as laser
        :return:
        """
        assert self.passQueue.empty()
        assert self.trajQueue.empty()


        self.laserObject = obj
        self.laserObject.generate_trajectory(self.startTime, self.endTime, self.steps)

        ## Define laser pass variables
        self._laserPasses = pd.DataFrame([], columns=self.dfCols)

    def gen_trajectories(self):
        """
        For all objects in queue, generate trajectory
        :return:
        """
        while True:
            try:
                obj = self.trajQueue.get()
                obj.generate_trajectory(self.startTime, self.endTime, self.steps)
                #self.message('Trajectory generated for {}'.format(obj.satName))
            except PropagationError as e:
                self.message('Got Propagation Error: {}'.format(e.msg))
                self.badTrajs.append(obj)
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
                relPosition = obj.trajectoryPos - self.laserObject.trajectoryPos
                relVelocity = obj.trajectoryVeloc - self.laserObject.trajectoryVeloc

                relDist = np.linalg.norm(relPosition, axis=1)
                indices = np.where(relDist < self.laserRange)
                smallDist = relDist[indices]/1000 #convert to km

                if len(smallDist) > 0:
                    tempDf = pd.DataFrame(columns=self.dfCols)
                    tempDf['Distance [km]'] = smallDist
                    tempDf['Time'] = [obj.trajectoryTimes[i] for i in indices[0]]
                    tempDf['Object'] = obj.satName
                    tempDf['Number'] = obj.satNum

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
    args = vars(parser.parse_args())

    parse = args['parse_data']
    dataDir = 'Data'
    savefile = 'satcat3_LEO.pickle'
    if parse or savefile not in os.listdir(os.path.join(os.getcwd(),dataDir)):
        print('Parsing tle file...')
        objects = parse_catalog(3)

        ## Throw out non-LEO objects
        nonLEO = []
        for i in range(0, len(objects)):
            obj = objects[i]
            h = (obj.a - Re) / 1e3  # km

            if h > 2000 or not obj.deb:
                nonLEO.append(i)

        for i in sorted(nonLEO, reverse=True):
            del objects[i]

        with open(os.path.join(os.getcwd(),dataDir,savefile), 'wb') as f:
            pickle.dump(objects, f)
    else:
        print('Loading pickle file')
        with open(os.path.join(os.getcwd(),dataDir,savefile), 'rb') as f:
            objects = pickle.load(f)

    ## Check that folders are made
    if not 'Plots' in os.listdir(os.getcwd()):
        os.makedirs(os.getcwd() + '/Plots')

    if not 'Data' in os.listdir(os.getcwd()):
        os.makedirs(os.getcwd() + '/Data')

    ######################
    ### Run Simulation ###
    ######################

    print('Starting sim...')
    print('Running with {} pieces of debris'.format(len(objects)))
    numDays = 1
    startTime = dt.datetime(2019,3,27,17,00,00)
    endTime = startTime + dt.timedelta(days=numDays)
    steps = numDays*24

    simulator = Simulator(startTime, endTime, steps)

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

    with open(os.path.join(os.getcwd(), dataDir, 'laserTLEObjects.pickle'), 'rb') as f:
        laser_objs = pickle.load(f)

    print('Running from {} to {}, with {} steps'.format(
        startTime.strftime(frmat), endTime.strftime(frmat), steps
    ))
    for laser_obj in laser_objs:
        simulator.setLaser(laser_obj)

        for obj in objects:
            if obj not in simulator.badTrajs:
                simulator.passQueue.put(obj)

        simulator.passQueue.join()

        ###############################
        ### Output Results and Plot ###
        ###############################

        passes = simulator.get_passes()
        print('')
        print('Laser Object: {}'.format(laser_obj.satName))
        uniquePasses = set([passes.iloc[i]['Number'] for i in range(0,len(passes))])
        print('Got {} passes, {} unique objects'.format(len(passes), len(uniquePasses)))
        for i in range(0, len(passes)):
            row = passes.iloc[i]
            print('{}: {} seen {} km away'.format(
                row['Time'].strftime(frmat),row['Object'], round(row['Distance [km]'],2)))

        output_file('Plots/laserPasses_{}.html'.format(simulator.laserObject.satNum))
        source = ColumnDataSource(passes)
        p = figure(title='Laser Passes Over Time', x_axis_label='Time', x_axis_type='datetime',y_axis_label='Distance [km]',
                   tooltips=[('Time','$x'),('Distance','$y'),('Object','@Object')])
        p.scatter(x='Time',y='Distance [km]',source=source)