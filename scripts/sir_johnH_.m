% ------------------------------------------------------------------------------
% 
% inversion of the SIR model on John Hopkins data
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
data_sir;
% ------------------------------------------------------------------------------
% initial conditions
corona_.nP = nP;                      % population
corona_.So = corona_.uo(1,1);         % susceptible
corona_.Io = corona_.uo(1,2);         % infectious
corona_.Ro = corona_.uo(1,3);         % removed
% ------------------------------------------------------------------------------
fprintf('\n\n   the data gives a SIR model:\n\n')
fprintf('       So susceptible   = %d \n' , corona_.So )
fprintf('       Io infected      = %d \n' , corona_.Io )
fprintf('       Ro removed       = %d \n' , corona_.Ro )
% ------------------------------------------------------------------------------
%%{
% ------------------------------------------------------------------------------
% 
%                            inversion
% 
% ------------------------------------------------------------------------------
b = 2.59e-01 ;
g = 2.10e-02 ;
%
p   = zeros(2,1);
p(1)= b;              % infection rate
p(2)= g;              % rate of removed
% ------------------------------------------------------------------------------
corona_.p = p;
corona_.p_ubounds = [2 2];
corona_.p_lbounds = zeros(1,size(p,1));
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = sir_fwd(corona_);
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
fprintf('       g removed rate                       = %2.2d \n' , g )
% ------------------------------------------------------------------------------
% 
%                     inversion routine
% 
% ------------------------------------------------------------------------------
corona_ = sir_inv_(corona_);
% ------------------------------------------------------------------------------
% run the fwd model again to see how we did
corona_.f = zeros(numel(corona_.t),1);
corona_.u = sir_fwd(corona_);
corona_.d = sir_u2d(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SIR model has:\n\n')
fprintf('       b infection rate                     = %2.2d \n' ,corona_.p(1) )
fprintf('       g removed rate                       = %2.2d \n' ,corona_.p(2) )
fprintf('\n\n\n')
%%{
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
corona_.u = sir_fwd(corona_);
corona_.d = sir_u2d(corona_);
t_obs = datetime(datestr(corona_.t_(:) + dates(istart)));
t     = datetime(datestr(corona_.t(:)  + dates(istart)));
% ------------------------------------------------------------------------------
figure;
hold on
plot(t, corona_.d(:,1),'k-')
plot(t, corona_.d(:,2),'r-')
plot(t_obs,corona_.do_(:,1),'ko','markersize',5)
plot(t_obs,corona_.do_(:,2),'ro','markersize',5)
hold off
axis tight
legend('Removed','Active')
xlabel('Time (days)')
ylabel('Population')
title(strcat('SIR model',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
%}