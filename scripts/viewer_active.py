import sys
sys.path.append('../graphics/')
from fancy_figure import fancy_figure
import numpy as np
import pandas as pd
# ------------------------------------------------------------------------------
guarda_path='../pics/'
path_      ='../data/tmp/'
# ------------------------------------------------------------------------------
dpi = 120
size=[9.6,5.76]
# ------------------------------------------------------------------------------
place_name=np.loadtxt(path_+'place_name.txt',dtype='str')
infec=np.loadtxt(path_+'time_series_covid19_confirmed_global.txt',dtype='float')
dead=np.loadtxt(path_+'time_series_covid19_deaths_global.txt',dtype='float')
reco=np.loadtxt(path_+'time_series_covid19_recovered_global.txt',dtype='float')
dates=np.loadtxt(path_+'dates.txt',dtype='str')
# ------------------------------------------------------------------------------
active = infec - ( reco+dead )
dates  = pd.to_datetime(dates)
dates  = np.array(dates,dtype=np.datetime64) 
# ------------------------------------------------------------------------------
if place_name.size > 1:
  place_name_=''
  for i_ in range(place_name.size):
    place_name_ = place_name_ + ' ' +place_name[i_]
  place_name = place_name_
else:
  place_name=place_name.item()
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
plt_=fancy_figure(
figsize=size,
x=dates,
curve=active,
colop='r',
symbol='.-',
margin=(0.03,0.03),
title='Active COVID-19 cases in '+place_name,
ylabel='Number of people',
xlabel='Time (days)',
fig_name=place_name+'-active',
guarda_path=guarda_path,
guarda=dpi).plotter_date()
# ------------------------------------------------------------------------------
plt_.show()
# ------------------------------------------------------------------------------







