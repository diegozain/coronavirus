function [istart,Io,Ro,Do] = start_date(I,R,D,dates,cut_off)

ndates = numel(dates);
istart= 1;
% while or(or(I(istart)==0 , R(istart)==0) , D(istart)==0)
while I(istart)<cut_off
 istart=istart+1;
end
Io = I(istart);
Ro = R(istart);
Do = D(istart);
end