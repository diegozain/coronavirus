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
% ------------------------------------------------------------------------------
% regularization term and weight
corona_.po= corona_.p;
corona_.bo= 0;
% ------------------------------------------------------------------------------
%                      inversion parameters
% ------------------------------------------------------------------------------
corona_.niter     = 1e+4;
corona_.tol_error = 1e-10;
% ------------------------------------------------------------------------------
% step-size
% 
% if you want a fixed step size,   uncomment corona_.step_ = ;
% if you want a parabolic search,  uncomment corona_.step_ = ;
% if you want an 1/norm(gradient), uncomment corona_.step_ave  = 1;
% ------------------------------------------------------------------------------
% corona_.step_ = 1e-10; % 1e-10
% corona_.perturb   = [-1e-5; 1e-3; 1e-1];
corona_.step_ave  = 1;
% ------------------------------------------------------------------------------
% bfgs stuff (not active unless uncommented in sir_inv.m)
% corona_.boH   = 1e-4;
corona_.H     = eye(numel(corona_.p));
corona_.p_    = zeros(size(corona_.p));
corona_.grad__= zeros(size(corona_.p));
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = sir_fwd(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n   our initial parameters are:\n\n')
fprintf('       b infection rate                     = %2.2d \n' , b )
fprintf('       g removed rate                       = %2.2d \n' , g )
% ------------------------------------------------------------------------------
% 
%                     inversion routine
% 
% ------------------------------------------------------------------------------
corona_ = sir_inv(corona_);
% ------------------------------------------------------------------------------
u = corona_.u;
p_history = corona_.p_history;
s_history = corona_.s_history;
E_history = corona_.E_history;
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SIR model has:\n\n')
fprintf('       b infection rate                     = %2.2d \n' ,corona_.p(1) )
fprintf('       g removed rate                       = %2.2d \n' ,corona_.p(2) )
fprintf('\n\n\n')
% ------------------------------------------------------------------------------
% get data d from field u
d = sir_u2d(corona_);
% ------------------------------------------------------------------------------
% plot history of the inversion
% ------------------------------------------------------------------------------
figure;
semilogy(p_history.','.-','markersize',20)
title('Parameter history')
xlabel('Iteration #')
ylabel('Value')
legend('b','g')
axis tight
simple_figure()

figure;
semilogy(E_history,'r.-','markersize',20)
title('Objecitve function history')
xlabel('Iteration #')
ylabel('Objecitve function')
axis tight
simple_figure()

figure;
semilogy(s_history,'b.-','markersize',20)
title('Step size history')
xlabel('Iteration #')
ylabel('Step size')
axis tight
simple_figure()
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
title('State vector, true (o), recovered (-)')
simple_figure();
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t, corona_.d(:,1),'k-')
plot(corona_.t, corona_.d(:,2),'r-')
plot(corona_.t_,corona_.do_(:,1),'ko','markersize',5)
plot(corona_.t_,corona_.do_(:,2),'ro','markersize',5)
hold off
axis tight
legend('Removed','Active')
xlabel('Time (days)')
ylabel('Population')
title(strcat('SIR model',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
%}