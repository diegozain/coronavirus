% ------------------------------------------------------------------------------
% 
% forward model of sir.
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
corona_ = struct;
% ------------------------------------------------------------------------------
% parameters
b = 0.9;
g = 0.1;
%
p   = zeros(2,1);
p(1)= b;           % infection rate
p(2)= g;           % rate of removed
corona_.p = p;
% ------------------------------------------------------------------------------
% time
corona_.dt= 1;
corona_.t = 0:corona_.dt:100;
% ------------------------------------------------------------------------------
% initial conditions
corona_.nP= 100;
corona_.Io = 2;
corona_.Ro = 0;
corona_.So = corona_.nP - (corona_.Io+corona_.Ro);
% ------------------------------------------------------------------------------
fprintf('\n\n   the true SEIR model has initial conditions:\n\n')
fprintf('       So susceptible   = %d \n' , corona_.So )
fprintf('       Io infected      = %d \n' , corona_.Io )
fprintf('       Ro recovered     = %d \n' , corona_.Ro )
fprintf('\n\n   with parameters:\n\n')
fprintf('       b infection rate  = %2.2d \n' , b )
fprintf('       g rate of removed = %2.2d \n' , g )
% ------------------------------------------------------------------------------
% 
%                               forward model
% 
% ------------------------------------------------------------------------------
corona_.f  = zeros(numel(corona_.t),1);
uo         = sir_fwd(corona_);
corona_.uo = uo;
% ------------------------------------------------------------------------------
% plot the true model
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t,uo(:,2),'r-','markersize',20)
plot(corona_.t,uo(:,3),'k-','markersize',20)
hold off
legend('infected','removed')
xlabel('Time')
ylabel('Population')
title('True')
simple_figure();
% ------------------------------------------------------------------------------
% 
% 
%                            inversion part
% 
% 
% ------------------------------------------------------------------------------
% initial guess & inversion parameters
b = 0.5 ; % 0.5 ; % 0.8 ;
g = 0.4 ; % 0.4 ; % 0.2 ;
%
p   = zeros(2,1);
p(1)= b;              % infection rate
p(2)= g;              % rate of removed
corona_.p = p;
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = sir_fwd(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the initial SIR model has:\n\n')
fprintf('       b infection rate  = %2.2d \n' , b )
fprintf('       g rate of removed = %2.2d \n' , g )
% ------------------------------------------------------------------------------
% see the adjoint
corona_.err = uo-u_init;
corona_.u   = u_init;
adj         = sir_adj(corona_);
% ------------------------------------------------------------------------------
% plot the adjoint
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t,adj(:,2),'r-','markersize',20)
plot(corona_.t,adj(:,3),'k-','markersize',20)
hold off
legend('infected','removed')
xlabel('Time')
ylabel('Population')
title('Adjoint')
simple_figure();
%%{
% ------------------------------------------------------------------------------
% 
%                          inversion routine
% 
% ------------------------------------------------------------------------------
%                        inversion parameters
% ------------------------------------------------------------------------------
corona_.niter     = 500; % 50000
corona_.tol_error = 1e-10;
% ------------------------------------------------------------------------------
% for step size
% 
% if you want a fixed step size, uncomment corona_.step_ = ;
% otherwise, if you want a parabolic search, comment corona_.step_ = ;
% 
corona_.step_     = 1e-4; % 3e-4
% corona_.perturb   = [-1e-8; 1e-8; 1e-7; 1e-6; 1e-4];
% corona_.perturb   = [-1e-5; 1e-3; 1e-2; 5e-2; 1e-2];
% ------------------------------------------------------------------------------
% inversion
corona_ = sir_inv(corona_);
% ------------------------------------------------------------------------------
u = corona_.u;
p_history = corona_.p_history;
s_history = corona_.s_history;
E_history = corona_.E_history;
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SEIR model has:\n\n')
fprintf('       b infection rate  = %2.2d \n' ,corona_.p(1) )
fprintf('       g rate of removed = %2.2d \n' ,corona_.p(2) )
% ------------------------------------------------------------------------------
% plot history of the inversion
% ------------------------------------------------------------------------------
figure;
plot(p_history.','.-','markersize',20)
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
% plot the true, initial and recovered models.
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t,uo(:,2),'r-','markersize',20)
plot(corona_.t,uo(:,3),'k-','markersize',20)
plot(corona_.t,u_init(:,2),'r-.','markersize',20)
plot(corona_.t,u_init(:,3),'k-.','markersize',20)
plot(corona_.t,u(:,2),'r.','markersize',20)
plot(corona_.t,u(:,3),'k.','markersize',20)
hold off
legend('infected','removed')
xlabel('Time')
ylabel('Population')
title('True ( - ) Initial ( -. ) Recovered ( . )')
simple_figure();
% ------------------------------------------------------------------------------
%}