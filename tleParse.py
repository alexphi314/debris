## Alex Philpott
## Loop through TLEs and print those that fall within certain altitude and inclination bounds
## In the end, not used to select a laser object, but useful for preliminary selection of test 'laser' satellites.

import argparse
import os
import pickle
import math

from astroUtils import parse_catalog, Re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Input arguments")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--parse_data", "-p", help="Re-parse sat cat file, rather than read stored data",
                          required=False, action="store_true", default=False)
    args = vars(parser.parse_args())

    parse = args['parse_data']
    dataDir = 'Data'
    savefile = 'satcat3_LEO_sat.pickle'
    if parse or savefile not in os.listdir(os.path.join(os.getcwd(),dataDir)):
        print('Parsing tle file...')
        objects = parse_catalog(3)

        ## Throw out non-LEO objects and all deb
        nonLEO = []
        for i in range(0, len(objects)):
            obj = objects[i]
            h = (obj.a - Re) / 1e3  # km

            if h > 2000 or obj.deb:
                nonLEO.append(i)

        for i in sorted(nonLEO, reverse=True):
            del objects[i]

        with open(os.path.join(os.getcwd(),dataDir,savefile), 'wb') as f:
            pickle.dump(objects, f)
    else:
        print('Loading pickle file')
        with open(os.path.join(os.getcwd(),dataDir,savefile), 'rb') as f:
            objects = pickle.load(f)

    target_hs = [700, 750, 800, 850, 900, 950, 1000]
    target_is = [0, 30, 45, 60, 75, 90, 96.6]
    h_band = 5
    i_band = 5

    for alt in target_hs:
        print('Target altitude: {} km'.format(alt))
        for inc in target_is:
            print('  Target inclination: {} deg'.format(inc))
            for obj in objects:
                h = (obj.a - Re)/1000
                i = math.degrees(obj.i)

                if h <= alt+h_band and h >=alt-h_band and i <= inc+i_band and i >= inc-i_band:
                    print('    {} ({}): alt {} km inc {} deg ecc {}'.format(
                        obj.satName, obj.satNum, round(h,1), round(i,1), round(obj.e,2)
                    ))

    satCatIDs = [309, 29050, 27560, 41181, 42806, 21576, 29047, 3510, 22219, 26536, 11971, 4256, 6126, 38744, 10358,
                 11080, 13492, 6909, 27001]
    satCatObjs = []
    for obj in objects:
        if obj.satNum in satCatIDs:
            satCatObjs.append(obj)

    with open(os.path.join(os.getcwd(), dataDir, 'laserTLEObjects.pickle'), 'wb') as f:
        pickle.dump(satCatObjs, f)