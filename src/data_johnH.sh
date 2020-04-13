#!/bin/bash
# ------------------------------------------------------------------------------
# get data global 
wget --no-check-certificate --content-disposition \
    https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv \
    -O ../data/time_series_covid19_confirmed_global.csv
#
wget --no-check-certificate \
    https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv \
    -O ../data/time_series_covid19_deaths_global.csv
#
wget --no-check-certificate \
    https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv \
    -O ../data/time_series_covid19_recovered_global.csv
# ------------------------------------------------------------------------------
# get data USA
wget --no-check-certificate --content-disposition \
    https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv \
    -O ../data/time_series_covid19_confirmed_US.csv
#
wget --no-check-certificate --content-disposition \
    https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv \
    -O ../data/time_series_covid19_deaths_US.csv
# ------------------------------------------------------------------------------