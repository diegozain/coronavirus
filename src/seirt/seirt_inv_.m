function corona_ = seirt_inv_(corona_)
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
Eo = corona_.Eo;
Io = corona_.Io;
Qo = corona_.Qo;
Ro = corona_.Ro;
Do = corona_.Do;
Po = corona_.Po;
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
fwd_ = @seirt_fwd_; 
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
    function d = seirt_fwd_(p,t)
    % --------------------------------------------------------------------------
    nt = numel(t);
    dt = mean(diff(t));
    % --------------------------------------------------------------------------
    % forcing term is zero... for now? this denotes people being born
    f = zeros(nt,1);
    % --------------------------------------------------------------------------
    % initial conditions
    u = zeros(nt,7);
    u(1,1) = So;
    u(1,2) = Eo;
    u(1,3) = Io;
    u(1,4) = Qo;
    u(1,5) = Ro;
    u(1,6) = Do;
    u(1,7) = Po;
    % --------------------------------------------------------------------------
    % Runge-Kutta of order 4
    % --------------------------------------------------------------------------
    for it=1:(nt-1);
     A = seirt_A_(u,p,t,nP,it);
     % runge-kutta
     k_1 = A*u(it,:).' + f(it,:).';
     k_2 = A*( u(it,:).' + 0.5*dt*k_1 ) + f(it,:).';
     k_3 = A*( u(it,:).' + 0.5*dt*k_2 ) + f(it,:).';
     k_4 = A*( u(it,:).' + dt*k_3 )     + f(it,:).';
     % output
     u(it+1,:) = u(it,:) + (1/6)*dt*(k_1 + 2*k_2 + 2*k_3 + k_4).';
    end
    % --------------------------------------------------------------------------
    d = seirt_u2d_(u,Mu2d);
    end
% ------------------------------------------------------------------------------
    function A = seirt_A_(u,p,t,nP,it)
    % --------------------------------------------------------------------------
    b = p(1);                % infection rate
    a = p(2);                % protection rate 
    g = p(3);                % inverse of average latent time
    d = p(4);                % inverse of average quarantine time
    l1= p(5);                % cure rate
    l2= p(6);                % cure rate
    k1= p(7);                % mortality
    k2= p(8);                % mortality 
    % --------------------------------------------------------------------------
    l = l1 * (1-exp(-l2*t(it)));
    k = k1 * ( exp(-k2*t(it)) );
    % --------------------------------------------------------------------------
    A      = zeros(7,7);
    % --------------------------------------------------------------------------
    A(1,1) = -(b*u(it,3)/nP) - a;
    A(2,1) =  b*u(it,3)/nP;
    A(7,1) =  a;
    A(2,2) = -g;
    A(3,2) =  g;
    A(3,3) = -d;
    A(4,3) =  d;
    A(4,4) = -(l+k);
    A(5,4) =  l;
    A(6,4) =  k;
    % --------------------------------------------------------------------------
    end
% ------------------------------------------------------------------------------
    function d = seirt_u2d_(u,Mu2d)
    nt= size(u,1);
    d = zeros(nt,4);
    for it=1:nt
     d_ = Mu2d * (u(it,:)).' ;
     d(it,:) = d_.';
    end
    end
% ------------------------------------------------------------------------------
end