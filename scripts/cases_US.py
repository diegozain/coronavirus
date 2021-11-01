import sys
sys.path.append('../graphics/')
from fancy_figure import fancy_figure
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
# ------------------------------------------------------------------------------
# get state name
state_= input ("state name  (eg Alabama)   : ")
# ------------------------------------------------------------------------------
# set data path
path_ = '../data/'
# ------------------------------------------------------------------------------
names_=['time_series_covid19_confirmed_US.csv','time_series_covid19_deaths_US.csv']
# ------------------------------------------------------------------------------
for name_ in names_:
    name__=path_+name_
    # --------------------------------------------------------------------------
    # load data
    data = pd.read_csv(name__,sep=',')
    # --------------------------------------------------------------------------
    # get from state
    data = data.loc[data['Province_State'] == state_]
    # --------------------------------------------------------------------------
    if name_=='time_series_covid19_confirmed_US.csv':
        # slice the data up to get the actual values and remove all the metadata stuff
        #  last two elements are "out of state" and "Unassigned" and I do not know what that means, so I'll remove them. 
        data = data.iloc[0:-2,11:]
        # get dates as a list of strings
        dates= data.columns.values
        dates= pd.to_datetime(dates)
        dates= np.array(dates,dtype=np.datetime64) 
        # get values for all counties and then sum them up
        data = data.values
        data = np.sum(data,axis=0)
        infec= data
    # --------------------------------------------------------------------------
    elif name_=='time_series_covid19_deaths_US.csv':
        npop = data.iloc[0:-2,11]
        data = data.iloc[0:-2,12:]
        # get values for all counties and then sum them up
        data = data.values
        data = np.sum(data,axis=0)
        dead = data
        npop = npop.values
        npop = np.sum(npop,axis=0)
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
infec= np.convolve(infec,np.ones(2),'same')
dead = np.convolve(dead,np.ones(2),'same')
# ------------------------------------------------------------------------------
infec= np.diff(infec)
dead = np.diff(dead)
dates= dates[0:-1]
# ------------------------------------------------------------------------------
guarda_path = '../pics/cases/'
dpi = 120
size=[9.6,5.76]
# ------------------------------------------------------------------------------
plt_=fancy_figure(
figsize=size,
x=dates,
curve=dead,
colop='k',
symbol='.-',
margin=(0.03,0.03),
holdon='on').plotter_date()
# ------------------------------------------------------------------------------
plt_=fancy_figure(
figsize=size,
x=dates,
curve=infec,
colop='r',
symbol='.-',
margin=(0.03,0.03),
title='Daily COVID-19 cases in '+state_,
ylabel='Number of people',
xlabel='Time (days)',
legends=['Dead','Infected'],
# frameon=True,
# legend_coord=(time_[1],0.8*dead.max()),
# legend_coord=(time_[1],0.8*max(infectious.max(),dead.max())),
fig_name='US-'+state_,
guarda_path=guarda_path,
guarda=dpi).plotter_date()
# ------------------------------------------------------------------------------
plt_.show()
# ------------------------------------------------------------------------------
