import sys
sys.path.append('../graphics/')
from fancy_figure import fancy_figure
import numpy as np
# ------------------------------------------------------------------------------
midi_s  = 1.7
# ------------------------------------------------------------------------------
guarda_path='../pics/'
path_      ='../data/tmp/'
# ------------------------------------------------------------------------------
dpi = 120
size=[9.6,9.6*0.6]
# ------------------------------------------------------------------------------
place_name=np.loadtxt(path_+'place_name.txt',dtype='str')
place_name=place_name.item()
# ------------------------------------------------------------------------------
time_,_,_     =fancy_figure.bring(path_,'time_',0,0)
infectious,_,_=fancy_figure.bring(path_,'infectious',0,0)
dead,_,_      =fancy_figure.bring(path_,'dead',0,0)
it_infec,_,_  =fancy_figure.bring(path_,'it_infec',0,0)
it_infec      = it_infec[0,0]
# ------------------------------------------------------------------------------
# matlab2py bullshit
time_ = time_- 366
time_ = time_.squeeze(axis=1)
it_infec = it_infec-1
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
fancy_figure(
figsize=size,
x=time_,
curve=infectious,
colop='b',
symbol='-',
holdon='on').plotter_date()
# ------------------------------------------------------------------------------
fancy_figure(
x=time_,
curve=dead,
colop='k',
symbol='-',
holdon='on').plotter_date()
# ------------------------------------------------------------------------------
fancy_figure(
x=time_[it_infec],
colop=(0.1843,0.8784,0.2196),
holdon='on').vline()
# ------------------------------------------------------------------------------
fancy_figure(
xmin=time_[it_infec],
xmax=time_[-1],
colop=(0.9216,0.2510,0.2039),
title='Estimate for '+place_name,
ylabel='Number of people',
xlabel='Time (days)',
legends=['Infectious','Dead','5\% of peak infectious','Open?'],
legend_coord=(time_[1],1.1*dead.max()),
fig_name=place_name+'-end',
guarda_path=guarda_path,
guarda=dpi).vspan()
# ------------------------------------------------------------------------------






