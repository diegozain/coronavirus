% ------------------------------------------------------------------------------
% 
% inversion of the SEIRt model on John Hopkins data
% 
%
% ------------------------------------------------------------------------------
clear
close all
clc
% ------------------------------------------------------------------------------
pwd_ = pwd;
cd ../
dir_paths;
cd(pwd_);
% ------------------------------------------------------------------------------
% loads the data. outputs,
% do    = [R D I repmat(nP,[size(R,1),1])];
% R     =  recovered
% D     =  dead
% I     =  infected confirmed
% nP    =  total population
% dates =  in datenum format
% istart=  index where the cut-off is done. 
%          this ensures initial conditions are nnz.
data_part;
% ------------------------------------------------------------------------------
% builds measuring operators, 
% normalizes the data (optional) and 
% interpolates the data (optional)
data_seirt;
% ------------------------------------------------------------------------------
% initial conditions
corona_.nP = nP;                      % population
corona_.So = corona_.uo(1,1);         % susceptible
corona_.Eo = corona_.uo(1,2);         % exposed       (unkown)
corona_.Io = corona_.uo(1,3);         % infectious
corona_.Qo = corona_.uo(1,4);         % quarantined   (unkown)
corona_.Ro = corona_.uo(1,5);         % recovered
corona_.Do = corona_.uo(1,6);         % dead 
corona_.Po = corona_.uo(1,7);         % insusceptible (unkown)
% ------------------------------------------------------------------------------
fprintf('\n\n   the data gives a SEIRt model:\n\n')
fprintf('       So susceptible   = %d \n' , corona_.So )
fprintf('       Eo exposed       = %d \n' , corona_.Eo )
fprintf('       Io infected      = %d \n' , corona_.Io )
fprintf('       Qo quarantined   = %d \n' , corona_.Qo )
fprintf('       Ro recovered     = %d \n' , corona_.Ro )
fprintf('       Do dead          = %d \n' , corona_.Do )
fprintf('       Po insusceptible = %d \n' , corona_.Po )
% ------------------------------------------------------------------------------
%%{
% ------------------------------------------------------------------------------
% 
%                            inversion
% 
% ------------------------------------------------------------------------------
b = 1e+00 ; % 1e+00 ; % 2.00e+00 ;
a = 6e-02 ; % 6e-02 ; % 2.84e-01 ;
g = 2e-01 ; % 2e-01 ; % 6.92e-02 ;
d = 1e-01 ; % 1e-01 ; % 2.13e-01 ;
l1= 1e-01 ; % 1e-01 ; % 1.70e-02 ;
l2= 5e-02 ; % 5e-02 ; % 2.00e+00 ;
k1= 1e-01 ; % 1e-01 ; % 1.99e-02 ;
k2= 5e-02 ; % 5e-02 ; % 3.89e-02 ;
%
p   = zeros(8,1);
p(1)= b;              % infection rate
p(2)= a;              % protection rate 
p(3)= g;              % inverse of average latent time
p(4)= d;              % inverse of average quarantine time
p(5)= l1;
p(6)= l2;             % l = l1 * (1-exp(-l2.*t)); % cure rate
p(7)= k1;             % k = k1 * ( exp(-k2*t(it)) ); % mortality
p(8)= k2; 
% ------------------------------------------------------------------------------
corona_.p = p;
corona_.p_ubounds = [1 2 1 1 1 2 1 2];
corona_.p_lbounds = zeros(1,size(p,1));
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),7);
u_init    = seirt_fwd(corona_);
% ------------------------------------------------------------------------------
%                      inversion parameters
% ------------------------------------------------------------------------------
corona_.niter    = 1e+4;
corona_.tol_x    = 1e-10;
corona_.tol_error= 1e-10;
corona_.fun_eval = 800;
corona_.algorithm='trust-region-reflective';
% corona_.algorithm='levenberg-marquardt';
% ------------------------------------------------------------------------------
fprintf('\n\n   our initial parameters are:\n\n')
fprintf('       b  infection rate                    = %2.2d \n' , b )
fprintf('       a  protection rate                   = %2.2d \n' , a )
fprintf('       g  inverse of average latent time    = %2.2d \n' , g )
fprintf('       d  inverse of average quarantine time= %2.2d \n' , d )
fprintf('       l1 cure rate                         = %2.2d \n' , l1)
fprintf('       l2 cure rate                         = %2.2d \n' , l2)
fprintf('       k1 mortality                         = %2.2d \n' , k1)
fprintf('       k2 mortality                         = %2.2d \n' , k2)
% ------------------------------------------------------------------------------
% 
%                     inversion routine
% 
% ------------------------------------------------------------------------------
corona_ = seirt_inv_(corona_);
% ------------------------------------------------------------------------------
% run the fwd model again to see how we did
corona_.f = zeros(numel(corona_.t),1);
corona_.u = seirt_fwd(corona_);
corona_.d = seirt_u2d(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SEIRt model has:\n\n')
fprintf('       b  infection rate                    = %2.2d \n' ,corona_.p(1))
fprintf('       a  protection rate                   = %2.2d \n' ,corona_.p(2))
fprintf('       g  inverse of average latent time    = %2.2d \n' ,corona_.p(3))
fprintf('       d  inverse of average quarantine time= %2.2d \n' ,corona_.p(4))
fprintf('       l1 cure rate                         = %2.2d \n' ,corona_.p(5) )
fprintf('       l2 cure rate                         = %2.2d \n' ,corona_.p(6) )
fprintf('       k1 mortality                         = %2.2d \n' ,corona_.p(7) )
fprintf('       k2 mortality                         = %2.2d \n' ,corona_.p(8) )
fprintf('\n\n\n')
%%{
% ------------------------------------------------------------------------------
% plot the true, initial and recovered models.
% ------------------------------------------------------------------------------
figure;
semilogy(corona_.t_,corona_.do_(:,3),'r.','markersize',25)
hold on
semilogy(corona_.t_,corona_.do_(:,1),'b.','markersize',25)
semilogy(corona_.t_,corona_.do_(:,2),'k.','markersize',25)
semilogy(corona_.t,u_init(:,4),'r-.')
semilogy(corona_.t,u_init(:,5),'b-.')
semilogy(corona_.t,u_init(:,6),'k-.')
semilogy(corona_.t, corona_.d(:,3),'r-')
semilogy(corona_.t, corona_.d(:,1),'b-')
semilogy(corona_.t, corona_.d(:,2),'k-')
hold off
axis tight
legend('Active','Recovered','Dead')
xlabel('Time (days)')
ylabel('Number of people')
title(strcat('True ( o ) Initial ( -. ) Recovered ( - )',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
% 
%                           now some projections
% 
% ------------------------------------------------------------------------------
days      = 3*31;
corona_.t = 0:corona_.dt:(corona_.t(end) + days);
corona_.nt= numel(corona_.t);
% run the fwd model again to see how we did
% ------------------------------------------------------------------------------
% here you can change the source and have new things happening at later times
corona_.f = zeros(numel(corona_.t),7);
% corona_.f(60,:) = [0,0,100,-100,0,0,0];
% ------------------------------------------------------------------------------
corona_.u = seirt_fwd(corona_);
corona_.d = seirt_u2d(corona_);
t_obs = datetime(datestr(corona_.t_(:) + dates(istart)));
t     = datetime(datestr(corona_.t(:)  + dates(istart)));
% ------------------------------------------------------------------------------
figure;
hold on
plot(t, corona_.d(:,3),'r-')
plot(t, corona_.d(:,1),'b-')
plot(t, corona_.d(:,2),'k-')
plot(t_obs,corona_.do_(:,3),'r.','markersize',25)
plot(t_obs,corona_.do_(:,1),'b.','markersize',25)
plot(t_obs,corona_.do_(:,2),'k.','markersize',25)
hold off
axis tight
grid on
legend('Active','Recovered','Dead')
xlabel('Time (days)')
ylabel('Number of people')
title(strcat('SEIRt model',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
name_pic = strcat('../pics/',place_name,'-fit');
fig = gcf;
print(gcf,name_pic,'-dpng','-r20')
print(gcf,name_pic,'-dpng')
% ------------------------------------------------------------------------------
% figure;
% hold on
% plot(t, corona_.u(:,3),'c-','linewidth',5)
% plot(t, corona_.d(:,2),'k-','linewidth',5)
% hold off
% grid on
% axis tight
% legend('Infectious','Dead')
% xlabel('Time (days)')
% ylabel('Number of people')
% title(strcat('SEIRt model',{' '},place_name))
% simple_figure();
% ------------------------------------------------------------------------------
% now let's see the rate of change of infectious cases.
% 
% taking the infectious cases equal to people in the hospital,
% the rate of change of infectious cases equals new cases.
% 
% the assumption that 
% 
% infectious cases == people in the hospital
% 
% is taken from the similarity between infectious cases (shown here) and
% people in the hospital taken from the website:
% 
% http://covid19.healthdata.org/
% 
dti = dtu(corona_.u(:,3),corona_.dt);
% ------------------------------------------------------------------------------
% allegedely, until 2 weeks have passed without any new cases things 
% will be operational again
% ------------------------------------------------------------------------------
% find moment in time when no new cases are reported
% the term fix(15/corona_.dt) is there because sometimes there is an artifact 
% in corona_.u(:,3) at early times.
it_zero = find_zero(dti( fix(15/corona_.dt):end ) , 0);
it_zero = it_zero + fix(15/corona_.dt);
if it_zero>numel(corona_.t)
 it_zero=numel(corona_.t);
end
% ------------------------------------------------------------------------------
% find four weeks (2*14 days) after such a time has happened
t_open = corona_.t(it_zero) + 2*14;
t_open = datetime(datestr( t_open  + dates(istart) ));
% ------------------------------------------------------------------------------
% figure;
% hold on
% plot(t, dti,'-')
% % plot(t, dtu(dti,corona_.dt),'.-')
% plot([t(it_zero),t(it_zero)] , [min(dti),max(dti)] ,'--','linewidth',5)
% plot([t_open,t_open] , [min(dti),max(dti)] ,'--','linewidth',5)
% hold off
% grid on
% axis tight
% xlabel('Time (days)')
% ylabel('Number of people')
% title(strcat('SEIRt model',{' '},place_name,{' - new cases'}))
% simple_figure();
% ------------------------------------------------------------------------------
% but maybe it would be a better idea to open up until after there are only, 
% say, 1000 infectious.
[imax,imaxi] = max(corona_.u(:,3));
it_infec = find_zero( corona_.u((imaxi+0):end,3) , 1e-1*imax);
if it_infec==1
 it_infec=imaxi;
else
 it_infec = imaxi + it_infec - 1;
end
% ------------------------------------------------------------------------------
ymini = min([min(corona_.u(:,3)),min(corona_.d(:,2))]);
ymaxi = max([max(corona_.u(:,3)),max(corona_.d(:,2))]);
% ------------------------------------------------------------------------------
% figure('Renderer', 'painters', 'Position', [0 0 900 400]);
% hold on
% plot(t, corona_.u(:,3),'c-','linewidth',5)
% plot(t, corona_.d(:,2),'k-','linewidth',5)
% plot([t(it_zero),t(it_zero)] , [ymini,ymaxi] ,'--','linewidth',5)
% plot([t_open,t_open] , [ymini,ymaxi] ,'--','linewidth',5)
% plot([t(it_infec),t(it_infec)] , [ymini,ymaxi] ,'--','linewidth',5)
% hold off
% grid on
% axis tight
% legend('Infectious','Dead (cumulative)','Zero new cases','Month after zero new cases','10% of peak infectious left')
% xlabel('Time (days)')
% ylabel('Number of people')
% title(strcat('SEIRt model',{' '},place_name))
% simple_figure();
% ------------------------------------------------------------------------------
% this part is for plotting in python
dead       = corona_.d(:,2);
infectious = corona_.u(:,3);
quarantined= corona_.u(:,4);
time_      = datenum(t);
% ------------------------------------------------------------------------------
save(strcat(path_,'dead'),'dead')
save(strcat(path_,'infectious'),'infectious')
save(strcat(path_,'quarantined'),'quarantined')
save(strcat(path_,'time_'),'time_')
save(strcat(path_,'t'),'t')
save(strcat(path_,'it_infec'),'it_infec')
% ------------------------------------------------------------------------------
!python viewer_end.py
!python viewer_active.py
% ------------------------------------------------------------------------------
% print(gcf,name_pic,'-dpng','-r20')
%}