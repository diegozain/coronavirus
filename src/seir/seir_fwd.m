function u = seir_fwd(corona_)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% seir model forward:
% 
% b = p(1);                % infection rate
% a = p(2);                % protection rate 
% g = p(3);                % inverse of average latent time
% d = p(4);                % inverse of average quarantine time
% l = p(5);                % cure rate
% k = p(6);                % mortality 
% 
% 
% 
% 
% 
% ------------------------------------------------------------------------------
% parameters
p = corona_.p;
% ------------------------------------------------------------------------------
% time and Population #
t = corona_.t;
nP= corona_.nP;
% ------------------------------------------------------------------------------
% source term
f = corona_.f;
% ------------------------------------------------------------------------------
dt= corona_.dt;
nt= numel(t);
% ------------------------------------------------------------------------------
% initial conditions
u = zeros(nt,7);
u(1,1) = corona_.So;
u(1,2) = corona_.Eo;
u(1,3) = corona_.Io;
u(1,4) = corona_.Qo;
u(1,5) = corona_.Ro;
u(1,6) = corona_.Do;
u(1,7) = corona_.Po;
% ------------------------------------------------------------------------------
% Runge-Kutta of order 4
% ------------------------------------------------------------------------------
for it=1:(nt-1);
 A = seir_A(u,p,nP,it);
 % runge-kutta
 k_1 = A*u(it,:).' + f(it,:).';
 k_2 = A*( u(it,:).' + 0.5*dt*k_1 ) + f(it,:).';
 k_3 = A*( u(it,:).' + 0.5*dt*k_2 ) + f(it,:).';
 k_4 = A*( u(it,:).' + dt*k_3 )     + f(it,:).';
 % output
 u(it+1,:) = u(it,:) + (1/6)*dt*(k_1 + 2*k_2 + 2*k_3 + k_4).';
end
% ------------------------------------------------------------------------------
end