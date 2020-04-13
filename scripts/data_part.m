% ------------------------------------------------------------------------------
% % get data from the internet
% prompt = '\n\n    Do you want to download data from the internet (y or n):  ';
% down_d = input(prompt,'s');
% % ----------------------------------------------------------------------------
% if strcmp(down_d,'y')
%  !sh ../src/data_johnH.sh
% end
% % ----------------------------------------------------------------------------
% extract what you want
prompt = '\n\n    Do you want to process a different country than last time (y or n):  ';
down_d = input(prompt,'s');
% ----------------------------------------------------------------------------
if strcmp(down_d,'y')
 !python ../src/process_johnH.py
end
% ------------------------------------------------------------------------------
% load into matlab
path_ = '../data/tmp/';
% ------------------------------------------------------------------------------
names_={'time_series_covid19_confirmed_global.txt','time_series_covid19_deaths_global.txt','time_series_covid19_recovered_global.txt'};
% ------------------------------------------------------------------------------
I = dlmread(strcat(path_,names_{1})); % infected
D = dlmread(strcat(path_,names_{2})); % dead
R = dlmread(strcat(path_,names_{3})); % recovered
% ------------------------------------------------------------------------------
% dates are in a special format
dates_ = fopen(strcat(path_,'dates.txt'));
dates_ = textscan(dates_,'%s');
nt     = size(dates_{1},1);
dates  = zeros(nt,1);
for it=1:nt
 dates(it) = datenum(dates_{1}(it),'mm/dd/yy');
end
% ------------------------------------------------------------------------------
% total population
nP      = dlmread(strcat(path_,'npop.txt'));
% ------------------------------------------------------------------------------
% I = window_mean(I,2);
% D = window_mean(D,2);
% R = window_mean(R,2);
% ------------------------------------------------------------------------------
% find starting date of infected
% cut_off = 0.3*max(I);
cut_off = 40;
cut_off = round(0.1*max(I));
[istart,Io,Ro,Do] = start_date(I,R,D,dates,cut_off);
% ------------------------------------------------------------------------------
D = D(istart:end);
I = I(istart:end);
R = R(istart:end);
% ------------------------------------------------------------------------------
do= [R D I repmat(nP,[size(R,1),1])];
% ------------------------------------------------------------------------------
% get name of place for a cute plot
place_name = fileread(strcat(path_,'place_name.txt'));
% ------------------------------------------------------------------------------
