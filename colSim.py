## Run collision simulator

import math
import argparse
import pickle
import os
import datetime as dt

import numpy as np

from astroUtils import parse_catalog, MEW, Re


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
    objRad = 5 #m
    startTime = dt.datetime.now()
    duration = 3600 #s
    timeStep = 3*60

    print('Running with {} pieces of debris'.format(len(objects)))

    positions = np.zeros([len(objects),3])
    for deltaT in range(0, duration+timeStep, timeStep):
        time = startTime + dt.timedelta(seconds=deltaT)
        print('')
        print('Time is {}'.format(time.strftime('%Y-%m-%d %H:%M:%S')))

        ## Update Position
        for i in range(0, len(objects)):
            obj = objects[i]
            positions[i,:] = obj.get_eci_pos(time)

        ## Check for conjunction
        for i in range(0, len(objects)):
            obj = objects[i]
            position = positions[i,:]

            relPos = positions - position
            relDist = np.linalg.norm(relPos, axis=1)
            indices = np.where(relDist<5)
            conjunctionJunction = relDist[indices]

            if len(conjunctionJunction) > 1:
                print('COLLISION DETECTED')
                print('Between {}'.format(', '.join(objects[i].satName for i in indices[0])))

