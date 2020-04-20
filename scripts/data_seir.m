% ------------------------------------------------------------------------------
%            Transform data to a linear combination of the data
% ------------------------------------------------------------------------------
% this part is crucial. Even though the observed data is
% 
% do = [R D I nP], 
% 
% the data to be fitted will be,
% 
% do = Md2d * [R D I-R-D nP] = [recovered, dead, active, population]
% 
% this model assumes the 'active' are all quarantined, so active == quarantined.
% ------------------------------------------------------------------------------
Md2d      = zeros(4);
Md2d(1,:) = [ 1  0 0 0];
Md2d(2,:) = [ 0  1 0 0];
Md2d(3,:) = [-1 -1 1 0];
Md2d(4,:) = [ 0  0 0 1];
do        = (Md2d*do.').';
% ------------------------------------------------------------------------------
% % normalize data
% do = do/nP;
% nP = 1;
% ------------------------------------------------------------------------------
% build corona_ struct
% ------------------------------------------------------------------------------
corona_ = struct;
% ------------------------------------------------------------------------------
% time & data interpolation
% ------------------------------------------------------------------------------
T = numel(dates)-(istart-1);
T = T-1; % minus one because of indicies starting at 1 and not 0, and dt=1 day.
% observed time
dt_= 1;        % days
t_ = 0:dt_:T;
nt_= numel(t_);
% interped time
dt= 0.5;       % days
t = 0:dt:T;
nt= numel(t);
% ------------------------------------------------------------------------------
% data obs
corona_.do_= do;
% data interp
corona_.do = interp1(t_,do,t);
% ------------------------------------------------------------------------------
% observed time
corona_.nt_= nt_;
corona_.dt_= dt_; % days
corona_.t_ = t_;  % days
% interped time
corona_.nt= nt;
corona_.dt= dt; % days
corona_.t = t;  % days
% ------------------------------------------------------------------------------
% this matrix takes the sate-vector u (solved by our fwd model) back to the 
% data (the measuring operator),
% 
% d = Mu2d * u
% 
% u(:,1) = susceptible
% u(:,2) = exposed
% u(:,3) = infectious
% u(:,4) = quarantined
% u(:,5) = recovered
% u(:,6) = dead 
% u(:,7) = insusceptible
% 
% this matrix will be used in the inversion,
% 
% d = Mu2d * u
% e = d - do
% adjoint source = Mu2d.' * e
% 
% ------------------------------------------------------------------------------
Mu2d      = zeros(4,7);
Mu2d(1,5) = 1;
Mu2d(2,6) = 1;
Mu2d(3,4) = 1;
Mu2d(4,:) = ones(1,7);
% ------------------------------------------------------------------------------
corona_.Mu2d = Mu2d;
% ------------------------------------------------------------------------------
% build observed data uo from observed data do.
% 
% u(:,1) = susceptible
% u(:,2) = exposed
% u(:,3) = infectious
% u(:,4) = quarantined
% u(:,5) = recovered
% u(:,6) = dead 
% u(:,7) = insusceptible
% 
M  = zeros(7,4);
% ------------------------------------------------------------------------------
% the idea for M is to formalize that the data is really R, D, A and nP, and 
% that the state-vector of our PDE has a linear combinations of them. 
% Let d be the vector [R; D; A; nP] for a moment in time. Then the state-vector 
% for a moment in time is,
% 
%   u = M*d
% 
% for example, if 
% M(2,:) = [1 1 1 0];
% then u(:,4), the exposed, will be = recovered + dead + active = infected
% ------------------------------------------------------------------------------
% these 3 are given
M(5,:) = [1   0 0 0];
M(6,:) = [0   1 0 0];
M(3,:) = [1   1 1 0];
% these are made up
M(4,:) = [0   0 1 0];
M(2,:) = [1   1 1 0];
M(7,:) = [0   0 0 0];
% this one should be exactly this... if the others are correct
M(1,:) = [0   0 0 1] - sum(M(2:7,:),1);
% ------------------------------------------------------------------------------
corona_.M  = M;
% ------------------------------------------------------------------------------
% collect observed data 
% ------------------------------------------------------------------------------
% not interpolated
nt = corona_.nt;
corona_.nt = corona_.nt_;
corona_.d  = corona_.do_;
corona_.uo_= seir_d2u(corona_);
corona_.nt = nt;
% interpolated
corona_.d  = corona_.do;
corona_.uo = seir_d2u(corona_);
% ------------------------------------------------------------------------------
