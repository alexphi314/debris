import math
import argparse
import pickle
import os

import numpy as np
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool, Legend, LogColorMapper, ColorBar, LogTicker
from bokeh.layouts import gridplot
from bokeh.palettes import Blues9

from .astroUtils import parse_catalog

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Input arguments")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--parse_data", "-p", help="Re-parse sat cat file, rather than read stored data",
                          required=False, action="store_true", default=False)
    #optional = parser.add_argument_group('optional arguments')
    args = vars(parser.parse_args())

    parse = args['parse_data']
    savefile = 'satcat3.pickle'
    if parse or savefile not in os.listdir(os.getcwd()):
        print('Parsing tle file...')
        objects = parse_catalog(3)
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
    upper_alt = 1000
    alt_delim = 5
    lower_i = 0
    upper_i = 180
    i_delim = 1

    np_alt = np.linspace(lower_alt, upper_alt, int((upper_alt+alt_delim-lower_alt)/alt_delim))
    altitudes = dict.fromkeys(list(np_alt),0)
    deb_altitudes = dict.fromkeys(list(np_alt),0)
    np_inc = np.linspace(lower_i, upper_i, int((upper_i+i_delim-lower_i)/i_delim))
    inclinations = dict.fromkeys(list(np_inc),0) #Range from 0 to 180 deg inclin. in 5 deg steps
    deb_inclinations = dict.fromkeys(list(np_inc),0)

    z = np.zeros([len(np_inc),len(np_alt)])
    deb_z = np.zeros([len(np_inc),len(np_alt)])

    x_vals = []
    y_vals = []
    deb_x_vals = []
    deb_y_vals = []

    debCounter = 0
    for obj in objects:
        if obj.deb:
            debCounter += 1

        inc = math.degrees(obj.i)
        alt = (obj.a - Re)/1000 #km

        nearest_alt = int(alt//alt_delim)*alt_delim
        nearest_inc = int(inc//i_delim)*i_delim

        if nearest_alt >= lower_alt and nearest_alt <= upper_alt:
            altitudes[nearest_alt]+= 1

            if obj.deb:
                deb_altitudes[nearest_alt]+=1

        if nearest_inc >= lower_i and nearest_inc <= upper_i:
            inclinations[nearest_inc]+=1

            if obj.deb:
                deb_inclinations[nearest_inc]+=1

        ## We need both within range to generate the image array
        if nearest_alt < lower_alt or nearest_alt > upper_alt \
            or nearest_inc < lower_i or nearest_inc > upper_i:
            continue

        x_coord = np.where(np_alt == nearest_alt)[0][0]
        y_coord = np.where(np_inc == nearest_inc)[0][0]

        z[y_coord][x_coord]+= 1
        x_vals.append(x_coord)
        y_vals.append(y_coord)

        if obj.deb:
            deb_z[y_coord][x_coord]+= 1
            deb_x_vals.append(x_coord)
            deb_y_vals.append(y_coord)

    print('{} TLEs, with {} deb'.format(len(objects),debCounter))
    Blues9.reverse()
    cmap = LogColorMapper(palette=Blues9,low=1,high=z.max())
    output_file('Plots/tle_distro.html')

    ## All TLE plots
    p = figure(title='Distribution of TLEs', x_axis_label='Altitude (km)', y_axis_label='TLE Amount',
               tooltips=[("x","$x"),("y","$y")])
    p.vbar(x=list(altitudes.keys()), top=list(altitudes.values()), bottom=0, width=alt_delim)

    p2 = figure(title='Distribution of TLEs', x_axis_label='Inclination (deg)', y_axis_label='TLE Amount',
                tooltips=[("x","$x"),("y","$y")])
    p2.vbar(x=list(inclinations.keys()), top=list(inclinations.values()), bottom=0, width=i_delim)

    p3 = figure(title='Distribution of TLEs',x_axis_label='Altitude (km)',y_axis_label='Inclination (deg)',
                x_range=(lower_alt,upper_alt),y_range=(lower_i,upper_i),
                tooltips=[("x","$x"),("y","$y"),("value","@image")])
    p3.image(image=[z],x=lower_alt,y=lower_i,dw=upper_alt-lower_alt,dh=upper_i-lower_i,color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap,location=(0,0),ticker=LogTicker())
    p3.add_layout(color_bar,'right')

    p4 = figure(title='Distribution of TLEs',x_axis_label='Altitude (km)',y_axis_label='Inclination (deg)',
                tooltips=[("x","$x"),("y","$y")])
    p4.scatter(x=x_vals,y=y_vals,marker='circle')

    ## Debris plots
    p5 = figure(title='Distribution of Debris', x_axis_label='Altitude (km)', y_axis_label='Debris Amount',
                tooltips=[("x","$x"),("y","$y")])
    p5.vbar(x=list(deb_altitudes.keys()), top=list(deb_altitudes.values()), bottom=0, width=alt_delim)

    p6 = figure(title='Distribution of Debris', x_axis_label='Inclination (deg)', y_axis_label='Debris Amount',
                tooltips=[("x","$x"),("y","$y")])
    p6.vbar(x=list(deb_inclinations.keys()), top=list(deb_inclinations.values()), bottom=0, width=i_delim)

    p7 = figure(title='Distribution of Debris', x_axis_label='Altitude (km)', y_axis_label='Inclination (deg)',
                x_range=(lower_alt, upper_alt), y_range=(lower_i, upper_i),
                tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    p7.image(image=[z], x=lower_alt, y=lower_i, dw=upper_alt - lower_alt, dh=upper_i - lower_i, color_mapper=cmap)
    color_bar = ColorBar(color_mapper=cmap, location=(0, 0), ticker=LogTicker())
    p7.add_layout(color_bar, 'right')

    p8 = figure(title='Distribution of Debris', x_axis_label='Altitude (km)', y_axis_label='Inclination (deg)',
                tooltips=[("x", "$x"), ("y", "$y")])
    p8.scatter(x=deb_x_vals, y=deb_y_vals, marker='circle')

    grid = gridplot([[p, p2],[p3, p4], [p5, p6], [p7, p8]])
    show(grid)