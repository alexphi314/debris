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

from astroUtils import parse_catalog, Re, PropagationError, Object

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
    tempLaserTLE = "LASER\n{}\n{}".format(
        '1   900U 64063C   19085.91254635 +.00000212 +00000-0 +21859-3 0  9999',
        '2   900 090.1517 021.2926 0026700 330.0700 101.6591 13.73206334708891'
    )
    laserRange = 100e3 #m
    def __init__(self, startTime, endTime, steps):
        """
        Init simulator manager class
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

        ##TODO: Use passed-in TLE for laser object definition
        self.laserObject = Object(self.tempLaserTLE)
        self.laserObject.generate_trajectory(self.startTime, self.endTime, self.steps)

        ## Define laser pass variables
        self.laserPasses = pd.DataFrame([],columns=['Time','Object','Distance [km]'])
        self.passQueue = queue.Queue()
        self.aLock = threading.Lock()

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
                self.passQueue.put(obj)
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
                relPosition = obj.trajectoryPos - self.laserObject.trajectoryPos
                relVelocity = obj.trajectoryVeloc - self.laserObject.trajectoryVeloc

                relDist = np.linalg.norm(relPosition, axis=1)
                indices = np.where(relDist < self.laserRange)
                smallDist = relDist[indices]/1000 #convert to km

                if len(smallDist) > 0:
                    tempDf = pd.DataFrame(columns=['Time','Object','Distance'])
                    tempDf['Distance [km]'] = smallDist
                    tempDf['Time'] = [obj.trajectoryTimes[i[0]] for i in indices if i is not None]
                    tempDf['Object'] = obj.satName

                    with self.aLock:
                        self.laserPasses = append(self.laserPasses, tempDf)

            except Exception as e:
                e_type, e_obj, e_tb = sys.exc_info()
                self.message('Got {} on line {}: {}'.format(e_type, e_tb.tb_lineno, e))
                pass
            self.passQueue.task_done()

    def message(self, msg):
        """
        Print message with thread name
        :param msg: message to print
        :return:
        """
        print('{}: {}'.format(threading.current_thread().name, msg))

if __name__ == "__main__":
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

    print('Starting sim...')
    print('Running with {} pieces of debris'.format(len(objects)))
    objRad = 5 #m
    startTime = dt.datetime(2019,3,27,17,00,00)
    endTime = startTime + dt.timedelta(days=2)
    steps = 48

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
    simulator.passQueue.join()
    print('Finished calculating passes')
    passes = simulator.laserPasses
    passes.to_csv('Data/laser_passes.csv',index=False)

    ## Check for conjunction
    # for i in range(0, len(objects)):
    #     obj = objects[i]
    #     position = positions[i,:]
    #
    #     relPos = positions - position
    #     relDist = np.linalg.norm(relPos, axis=1)
    #     indices = np.where(relDist<5)
    #     conjunctionJunction = relDist[indices]
    #
    #     if len(conjunctionJunction) > 1:
    #         print('COLLISION DETECTED')
    #         print('Between {}'.format(', '.join(objects[i].satName for i in indices[0])))

