function A = seirt_A(u,p,t,nP,it)
% -----------------------------------------------------------------------------
% diego domenzain
% spring 2020 while at Colorado School of Mines
% ------------------------------------------------------------------------------
% inversion of the SEIR model using the adjoint method.
% this function computes the A matrix.
% ------------------------------------------------------------------------------
b = p(1);                % infection rate
a = p(2);                % protection rate 
g = p(3);                % inverse of average latent time
d = p(4);                % inverse of average quarantine time
l1= p(5);                % cure rate
l2= p(6);                % cure rate
k1= p(7);                % mortality
k2= p(8);                % mortality 
% ------------------------------------------------------------------------------
l = l1 * (1-exp(-l2*t(it)));
k = k1 * ( exp(-k2*t(it)) );
% ------------------------------------------------------------------------------
A      = zeros(7,7);
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------
end