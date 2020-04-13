% ------------------------------------------------------------------------------
% 
% inversion of the SEIR model on John Hopkins data
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
data_seir;
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
fprintf('\n\n   the data gives a SEIR model:\n\n')
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
l = 1e-01 ; % 1e-01 ; % 1.70e-02 ;
k = 1e-01 ; % 5e-02 ; % 2.00e+00 ;
%
p   = zeros(6,1);
p(1)= b;              % infection rate
p(2)= a;              % protection rate 
p(3)= g;              % inverse of average latent time
p(4)= d;              % inverse of average quarantine time
p(5)= l;              % cure rate
p(6)= k;              % mortality
% ------------------------------------------------------------------------------
corona_.p = p;
corona_.p_ubounds = [1 2 1 1 1 1];
corona_.p_lbounds = zeros(1,size(p,1));
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = seir_fwd(corona_);
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
fprintf('       b infection rate                     = %2.2d \n' , b )
fprintf('       a protection rate                    = %2.2d \n' , a )
fprintf('       g inverse of average latent time     = %2.2d \n' , g )
fprintf('       d inverse of average quarantine time = %2.2d \n' , d )
fprintf('       l1 cure rate                         = %2.2d \n' , l)
fprintf('       k1 mortality                         = %2.2d \n' , k)
% ------------------------------------------------------------------------------
% 
%                     inversion routine
% 
% ------------------------------------------------------------------------------
corona_ = seir_inv_(corona_);
% ------------------------------------------------------------------------------
% run the fwd model again to see how we did
corona_.f = zeros(numel(corona_.t),1);
corona_.u = seir_fwd(corona_);
corona_.d = seir_u2d(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SEIR model has:\n\n')
fprintf('       b infection rate                     = %2.2d \n' ,corona_.p(1) )
fprintf('       a protection rate                    = %2.2d \n' ,corona_.p(2) )
fprintf('       g inverse of average latent time     = %2.2d \n' ,corona_.p(3) )
fprintf('       d inverse of average quarantine time = %2.2d \n' ,corona_.p(4) )
fprintf('       l cure rate                          = %2.2d \n' ,corona_.p(5) )
fprintf('       k mortality                          = %2.2d \n' ,corona_.p(6) )
fprintf('\n\n\n')
%%{
% ------------------------------------------------------------------------------
% plot the true, initial and recovered models.
% ------------------------------------------------------------------------------
figure;
semilogy(corona_.t_,corona_.do_(:,3),'ro','markersize',5)
hold on
semilogy(corona_.t_,corona_.do_(:,1),'bo','markersize',5)
semilogy(corona_.t_,corona_.do_(:,2),'ko','markersize',5)
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
ylabel('Population')
title('True ( o ) Initial ( -. ) Recovered ( - )')
simple_figure();
% ------------------------------------------------------------------------------
figure;
semilogy(corona_.t_,corona_.uo_,'o','markersize',5)
hold on
semilogy(corona_.t, corona_.u,'-','markersize',20)
hold off
axis tight
xlabel('Time (days)')
ylabel('Population')
title('State-vector, true (o), recovered (-)')
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
corona_.f = zeros(numel(corona_.t),1);
corona_.u = seir_fwd(corona_);
corona_.d = seir_u2d(corona_);
t_obs = datetime(datestr(corona_.t_(:) + dates(istart)));
t     = datetime(datestr(corona_.t(:)  + dates(istart)));
% ------------------------------------------------------------------------------
figure;
hold on
plot(t, corona_.d(:,3),'r-')
plot(t, corona_.d(:,1),'b-')
plot(t, corona_.d(:,2),'k-')
plot(t_obs,corona_.do_(:,3),'ro','markersize',5)
plot(t_obs,corona_.do_(:,1),'bo','markersize',5)
plot(t_obs,corona_.do_(:,2),'ko','markersize',5)
hold off
axis tight
legend('Active','Recovered','Dead')
xlabel('Time (days)')
ylabel('Population')
title(strcat('SEIR model',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
%}