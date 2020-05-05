import os
import pandas as pd
import numpy as np
# ------------------------------------------------------------------------------
# get country name
country_  = input ("country name  (case sensitive)            : ")
province_ = input ("province name (if none then press enter)  : ")
# ------------------------------------------------------------------------------
# set data path
path_ = '../data/'
# ------------------------------------------------------------------------------
# set names to read. maybe fix this for something more clever like ls?
names_=['time_series_covid19_confirmed_global.csv','time_series_covid19_deaths_global.csv','time_series_covid19_recovered_global.csv']
# ------------------------------------------------------------------------------
for name_ in names_:
    name__=path_+name_
    # --------------------------------------------------------------------------
    # load data
    data = pd.read_csv(name__,sep=',')
    # --------------------------------------------------------------------------
    # get from country
    data = data.loc[data['Country/Region'] == country_]
    # --------------------------------------------------------------------------
    # get from province
    if province_ != '':
        data = data.loc[data['Province/State'] == province_]
    # --------------------------------------------------------------------------
    # slice the data up to get the actual values and remove all the metadata stuff
    data = data.iloc[0:,4:]
    # get dates as a list of strings
    dates= data.columns.values
    # get values for all time and all provinces and then sum them up
    data = data.values
    data = np.sum(data,axis=0)
    # --------------------------------------------------------------------------
    os.chdir('../data/tmp')
    os.system('touch '+name_[:-4]+'.txt')
    np.savetxt(name_[:-4]+'.txt',data)
    os.system('touch dates.txt')
    np.savetxt('dates.txt',dates,fmt='%s')
    os.chdir('../../src/')
# ------------------------------------------------------------------------------
# finaly, get population number for that country or province
name_ = 'time_series_covid19_population_global.csv'
name__=path_+name_
npop  = pd.read_csv(name__,sep=',')
# --------------------------------------------------------------------------
# get from country
npop = npop.loc[npop['Country/Region'] == country_]
# --------------------------------------------------------------------------
# get from province
if province_ != '':
    npop = npop.loc[npop['Province/State'] == province_]
    place_name = country_ + ', ' + province_
else:
    place_name = country_
place_name = [place_name]
# ------------------------------------------------------------------------------
# print(npop.shape)
npop = npop.iloc[0,2]
os.chdir('../data/tmp')
os.system('touch npop.txt')
np.savetxt('npop.txt',[npop])
os.system('touch place_name.txt')
np.savetxt('place_name.txt',place_name,fmt="%s")
os.chdir('../../src/')
# ------------------------------------------------------------------------------
