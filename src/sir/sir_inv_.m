function corona_ = sir_inv_(corona_)
% ------------------------------------------------------------------------------
% inspired by https://github.com/ECheynet/SEIR
% 
% diego domenzain, spring 2020
% @ Colorado School of Mines
% ------------------------------------------------------------------------------
do        = corona_.do;       % observed data
nP        = corona_.nP;       % total population
Mu2d      = corona_.Mu2d;     % state-vector to data matrix (measuring operator)
p         = corona_.p;        % initial parameters
p_ubounds = corona_.p_ubounds;% upper bounds for p
p_lbounds = corona_.p_lbounds;% lower bounds for p
t         = corona_.t;        % time vector
% ------------------------------------------------------------------------------
% initial conditions
So = corona_.So;
Io = corona_.Io;
Ro = corona_.Ro;
% ------------------------------------------------------------------------------
tol_x    = corona_.tol_x;
tol_error= corona_.tol_error;
fun_eval = corona_.fun_eval;
algorithm= corona_.algorithm;
% ------------------------------------------------------------------------------
% 
%                   invoke the magic of lsqcurvefit
% 
% ------------------------------------------------------------------------------
% options for lsqcurvfit
options = optimset('TolX',tol_x,'TolFun',tol_error,'MaxFunEvals',fun_eval,'Display','iter','Algorithm',algorithm);
% ------------------------------------------------------------------------------
% transform the forward model (a nested function) into an anonymous function
fwd_ = @sir_fwd_; 
% ------------------------------------------------------------------------------
[p,~,e_,~,output,lambda,J_] = lsqcurvefit(@(p,t) fwd_(p,t),p,t,do,p_lbounds,p_ubounds,options);
% ------------------------------------------------------------------------------
corona_.p  = p;
corona_.e_ = e_;
corona_.J_ = J_;
corona_.output = output;
corona_.lambda = lambda;
% ------------------------------------------------------------------------------
% 
%                           nested fwd model
% 
% ------------------------------------------------------------------------------
    function d = sir_fwd_(p,t)
    % --------------------------------------------------------------------------
    nt = numel(t);
    dt = mean(diff(t));
    % --------------------------------------------------------------------------
    % forcing term is zero... for now? this denotes people being born
    f = zeros(nt,1);
    % --------------------------------------------------------------------------
    % initial conditions
    u = zeros(nt,3);
    u(1,1) = So;
    u(1,2) = Io;
    u(1,3) = Ro;
    % --------------------------------------------------------------------------
    % Runge-Kutta of order 4
    % --------------------------------------------------------------------------
    for it=1:(nt-1);
     A = sir_A_(u,p,t,nP,it);
     % runge-kutta
     k_1 = A*u(it,:).' + f(it,:).';
     k_2 = A*( u(it,:).' + 0.5*dt*k_1 ) + f(it,:).';
     k_3 = A*( u(it,:).' + 0.5*dt*k_2 ) + f(it,:).';
     k_4 = A*( u(it,:).' + dt*k_3 )     + f(it,:).';
     % output
     u(it+1,:) = u(it,:) + (1/6)*dt*(k_1 + 2*k_2 + 2*k_3 + k_4).';
    end
    % --------------------------------------------------------------------------
    d = sir_u2d_(u,Mu2d);
    end
% ------------------------------------------------------------------------------
    function A = sir_A_(u,p,t,nP,it)
    % --------------------------------------------------------------------------
    b = p(1);                % infection rate
    g = p(2);                % inverse of average latent time 
    % --------------------------------------------------------------------------
    A      = zeros(3,3);
    % --------------------------------------------------------------------------
    A(1,1) = -(b*u(it,2)/nP);
    A(2,1) =   b*u(it,2)/nP;
    A(2,2) = -g;
    A(3,2) =  g;
    % --------------------------------------------------------------------------
    end
% ------------------------------------------------------------------------------
    function d = sir_u2d_(u,Mu2d)
    nt= size(u,1);
    d = zeros(nt,3);
    for it=1:nt
     d_ = Mu2d * (u(it,:)).' ;
     d(it,:) = d_.';
    end
    end
% ------------------------------------------------------------------------------
end