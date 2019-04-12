import datetime as dt

from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.palettes import Blues9
from bokeh.layouts import gridplot

import pandas as pd

x_data = [.25,.5,3,4,5]
y_data = [1,2,3,4,5]

p = figure(title='Debris Decays', x_axis_label='Time',
                y_axis_label='Number of Debris Decays')
p.line(x=x_data,y=y_data)

#grid = gridplot([[p]])
show(p)