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
figure;
hold on
plot(corona_.t_,corona_.do_(:,3),'r-','markersize',20)
plot(corona_.t_,corona_.do_(:,1),'b-','markersize',20)
plot(corona_.t_,corona_.do_(:,2),'k-','markersize',20)
hold off
axis tight
legend('Active','Recovered','Dead')
xlabel('Time (days)')
ylabel('Population')
title('True')
simple_figure();
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
% corona_.step_ = 1e+5;
corona_.perturb   = [-1e-5; 1e-1; 1e+1];
% corona_.step_ave  = 1;
% ------------------------------------------------------------------------------
% bfgs stuff (not active unless uncommented in seirt_inv.m)
% corona_.boH   = 1e-4;
corona_.H     = eye(numel(corona_.p));
corona_.p_    = zeros(size(corona_.p));
corona_.grad__= zeros(size(corona_.p));
% ------------------------------------------------------------------------------
% record initial guess to see how far off we were
corona_.f = zeros(numel(corona_.t),1);
u_init    = seirt_fwd(corona_);
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
corona_ = seirt_inv(corona_);
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
% ------------------------------------------------------------------------------
%                   plot history of the inversion
% ------------------------------------------------------------------------------
p_history = corona_.p_history;
s_history = corona_.s_history;
E_history = corona_.E_history;
% ------------------------------------------------------------------------------
figure;
semilogy(p_history.','.-','markersize',20)
title('Parameter history')
xlabel('Iteration #')
ylabel('Value')
legend('b','a','g','d','l1','l2','k1','k2')
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
hold on
plot(corona_.t, corona_.d(:,3),'r-')
plot(corona_.t, corona_.d(:,1),'b-')
plot(corona_.t, corona_.d(:,2),'k-')
plot(corona_.t_,corona_.do_(:,3),'ro','markersize',5)
plot(corona_.t_,corona_.do_(:,1),'bo','markersize',5)
plot(corona_.t_,corona_.do_(:,2),'ko','markersize',5)
hold off
axis tight
legend('Active','Recovered','Dead')
xlabel('Time (days)')
ylabel('Population')
title(strcat('SEIRt model',{' '},place_name))
simple_figure();
% ------------------------------------------------------------------------------
%}