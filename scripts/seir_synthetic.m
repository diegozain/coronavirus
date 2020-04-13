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
a = 0.4;
g = 0.1;
d = 0.3;
l = 0.2;
k = 0.1;
%
p   = zeros(6,1);
p(1)= b;           % infection rate
p(2)= a;           % protection rate 
p(3)= g;           % inverse of average latent time
p(4)= d;           % inverse of average quarantine time
p(5)= l;           % cure rate
p(6)= k;           % mortality 
corona_.p = p;
% ------------------------------------------------------------------------------
% time
corona_.dt= 1;
corona_.t = 0:corona_.dt:100;
% ------------------------------------------------------------------------------
% initial conditions
corona_.nP= 100;
corona_.Eo = 0; % 1
corona_.Io = 2;
corona_.Qo = 0; % 3
corona_.Ro = 0;
corona_.Do = 0;
corona_.Po = 0;
corona_.So = corona_.nP - (corona_.Eo+corona_.Io+corona_.Qo+corona_.Ro+corona_.Do+corona_.Po);
% ------------------------------------------------------------------------------
fprintf('\n\n   the true SEIR model has initial conditions:\n\n')
fprintf('       So susceptible   = %d \n' , corona_.So )
fprintf('       Eo exposed       = %d \n' , corona_.Eo )
fprintf('       Io infected      = %d \n' , corona_.Io )
fprintf('       Qo quarantined   = %d \n' , corona_.Qo )
fprintf('       Ro recovered     = %d \n' , corona_.Ro )
fprintf('       Do dead          = %d \n' , corona_.Do )
fprintf('       Po insusceptible = %d \n' , corona_.Po)
fprintf('\n\n   with parameters:\n\n')
fprintf('       b infection rate                     = %2.2d \n' , b )
fprintf('       a protection rate                    = %2.2d \n' , a )
fprintf('       g inverse of average latent time     = %2.2d \n' , g )
fprintf('       d inverse of average quarantine time = %2.2d \n' , d )
fprintf('       l cure rate                          = %2.2d \n' , l )
fprintf('       k mortality                          = %2.2d \n' , k )
% ------------------------------------------------------------------------------
% 
%                               forward model
% 
% ------------------------------------------------------------------------------
corona_.f  = zeros(numel(corona_.t),1);
uo         = seir_fwd(corona_);
corona_.uo = uo;
% ------------------------------------------------------------------------------
% plot the true model
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t,uo(:,3),'r-','markersize',20)
plot(corona_.t,uo(:,5),'b-','markersize',20)
plot(corona_.t,uo(:,6),'k-','markersize',20)
hold off
legend('infected','recovered','dead')
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
a = 0.1 ; % 0.1 ; % 0.3 ;
g = 0.4 ; % 0.4 ; % 0.2 ;
d = 0.6 ; % 0.6 ; % 0.4 ;
l = 0.4 ; % 0.4 ; % 0.1 ;
k = 0.7 ; % 0.7 ; % 0.2 ;
%
p   = zeros(6,1);
p(1)= b;              % infection rate
p(2)= a;              % protection rate 
p(3)= g;              % inverse of average latent time
p(4)= d;              % inverse of average quarantine time
p(5)= l;              % cure rate
p(6)= k;              % mortality 
corona_.p = p;
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = seir_fwd(corona_);
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the initial SIR model has:\n\n')
fprintf('       b infection rate                     = %2.2d \n' , b )
fprintf('       a protection rate                    = %2.2d \n' , a )
fprintf('       g inverse of average latent time     = %2.2d \n' , g )
fprintf('       d inverse of average quarantine time = %2.2d \n' , d )
fprintf('       l cure rate                          = %2.2d \n' , l )
fprintf('       k mortality                          = %2.2d \n' , k )
% ------------------------------------------------------------------------------
% see the adjoint
corona_.err = uo-u_init;
corona_.u   = u_init;
adj         = seir_adj(corona_);
% ------------------------------------------------------------------------------
% plot the adjoint
% ------------------------------------------------------------------------------
figure;
hold on
plot(corona_.t,adj(:,3),'r-','markersize',20)
plot(corona_.t,adj(:,5),'b-','markersize',20)
plot(corona_.t,adj(:,6),'k-','markersize',20)
hold off
legend('infected','recovered','dead')
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
corona_.niter     = 5000; % 50000
corona_.tol_error = 1e-10;
% ------------------------------------------------------------------------------
% for step size
% 
% if you want a fixed step size, uncomment corona_.step_ = ;
% otherwise, if you want a parabolic search, comment corona_.step_ = ;
% 
corona_.step_     = 3e-4; % 3e-4
% corona_.perturb   = [-1e-8; 1e-8; 1e-7; 1e-6; 1e-4];
% corona_.perturb   = [-1e-5; 1e-3; 1e-2; 5e-2; 1e-2];
% corona_.step_ave  = 1;
% ------------------------------------------------------------------------------
% corona_.H     = eye(6);
% corona_.boH   = 1e+3;
% corona_.p_    = zeros(size(corona_.p));
% corona_.grad__= zeros(size(corona_.p));
% ------------------------------------------------------------------------------
% inversion
corona_ = seir_inv(corona_);
% ------------------------------------------------------------------------------
u = corona_.u;
p_history = corona_.p_history;
s_history = corona_.s_history;
E_history = corona_.E_history;
% ------------------------------------------------------------------------------
fprintf('\n\n ------------')
fprintf('\n\n   the recovered SEIR model has:\n\n')
fprintf('       b infection rate                     = %2.2d \n' ,corona_.p(1) )
fprintf('       a protection rate                    = %2.2d \n' ,corona_.p(2) )
fprintf('       g inverse of average latent time     = %2.2d \n' ,corona_.p(3) )
fprintf('       d inverse of average quarantine time = %2.2d \n' ,corona_.p(4) )
fprintf('       l cure rate                          = %2.2d \n' ,corona_.p(5) )
fprintf('       k mortality                          = %2.2d \n' ,corona_.p(6) )
% ------------------------------------------------------------------------------
% plot history of the inversion
% ------------------------------------------------------------------------------
figure;
plot(p_history.','.-','markersize',20)
title('Parameter history')
xlabel('Iteration #')
ylabel('Value')
legend('b','a','g','d','l','k')
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
plot(corona_.t,uo(:,3),'r-','markersize',20)
plot(corona_.t,uo(:,5),'b-','markersize',20)
plot(corona_.t,uo(:,6),'k-','markersize',20)
plot(corona_.t,u_init(:,3),'r-.','markersize',20)
plot(corona_.t,u_init(:,5),'b-.','markersize',20)
plot(corona_.t,u_init(:,6),'k-.','markersize',20)
plot(corona_.t,u(:,3),'r.','markersize',20)
plot(corona_.t,u(:,5),'b.','markersize',20)
plot(corona_.t,u(:,6),'k.','markersize',20)
hold off
legend('infected','recovered','dead')
xlabel('Time')
ylabel('Population')
title('True ( - ) Initial ( -. ) Recovered ( . )')
simple_figure();
% ------------------------------------------------------------------------------
%}