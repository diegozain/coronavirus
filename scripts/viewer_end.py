import sys
sys.path.append('../graphics/')
from fancy_figure import fancy_figure
import numpy as np
# ------------------------------------------------------------------------------
guarda_path='../pics/'
path_      ='../data/tmp/'
# ------------------------------------------------------------------------------
dpi = 120
size=[9.6,5.76]
# ------------------------------------------------------------------------------
place_name=np.loadtxt(path_+'place_name.txt',dtype='str')
# ------------------------------------------------------------------------------
if place_name.size > 1:
  place_name_=''
  for i_ in range(place_name.size):
    place_name_ = place_name_ + ' ' +place_name[i_]
  place_name = place_name_
else:
  place_name=place_name.item()
# ------------------------------------------------------------------------------
time_,_,_      =fancy_figure.bring(path_,'time_',0,0)
infectious,_,_ =fancy_figure.bring(path_,'infectious',0,0)
quarantined,_,_=fancy_figure.bring(path_,'quarantined',0,0)
dead,_,_       =fancy_figure.bring(path_,'dead',0,0)
it_infec,_,_   =fancy_figure.bring(path_,'it_infec',0,0)
it_infec       = it_infec[0,0]
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
x=time_,
curve=quarantined,
colop='m',
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
legends=['Infectious','Dead','Quarantined','10\% of peak infectious','Open?'],
frameon=True,
# legend_coord=(time_[1],0.8*dead.max()),
# legend_coord=(time_[1],0.8*max(infectious.max(),dead.max())),
fig_name=place_name+'-end',
guarda_path=guarda_path,
guarda=dpi).vspan()
# ------------------------------------------------------------------------------







