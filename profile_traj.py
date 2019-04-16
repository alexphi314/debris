## Profile code to speed up execution
import cProfile
import argparse
import os
import sys
import pickle
import threading
import datetime as dt
import pstats

from astroUtils import parse_catalog, Re, Laser

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Input arguments")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--parse_data", "-p", help="Re-parse sat cat file, rather than read stored data",
                          required=False, action="store_true", default=False)
    args = vars(parser.parse_args())

    print('Running with call: {}'.format(sys.argv[0:]))

    parse = args['parse_data']
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

    objects = deb + sats
    obj = objects[0]
    laser = Laser(sats[10],True)

    assert obj != laser

    matlabLock = threading.RLock()
    laserLock = threading.RLock()
    startTime = dt.datetime(2019,4,4,17,00,00)
    endTime = dt.datetime(2019,4,11,17,00,00)
    steps = 7*24

    laser.generate_trajectory(startTime,endTime,steps,matlabLock,laserLock)

    p = cProfile.Profile()
    p.run('obj.generate_trajectory(startTime,endTime,steps,matlabLock,laserLock,laser)')
    stats = pstats.Stats(p)
    stats.sort_stats('cumulative')
    stats.print_stats('astroUtils')

    for indx,obj in enumerate(objects):
        if obj.satNum == 43689:
            print('Laser period: {} min'.format(obj.T/60))


