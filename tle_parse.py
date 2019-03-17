import math
import datetime as dt

import bokeh

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

        ## Define times
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


if __name__ == "__main__":
    objects = parse_catalog()

