function [Rq,r] = esterr_kernel_lcp(q,p,GammaK,s,normV,tol,zbar,qbar,fellZA,z,dz)
% linear combination
n = q-2;
% Gautschi algorithm for kernel evaluations
nu0 = max(ceil(n+log(1/eps)./(2*log(abs((zbar+sqrt(zbar-1).*sqrt(zbar+1)))))));
ab = r_jacobi(nu0+1,1,1);
b = ab(:,2);
y0 = ones(size(zbar));
ym1 = zeros(size(zbar));
y = 1; % to handle q=2 (n=0, next for loop skipped)
for k = 0:n-1
  y = zbar.*y0-b(k+1)*ym1;
  ym1 = y0;
  y0 = y;
end
erre(nu0+2,:) = zeros(size(zbar));
for k = nu0:-1:0
  erre(k+1,:) = b(k+1)./(zbar-erre(k+2,:));
end
rho(1,:) = ones(size(zbar));
for k = 0:n
  rho(k+2,:) = erre(k+1,:).*rho(k+1,:);
end
Kn_Gautschi = rho(n+2,:)./y./(1-zbar.*zbar);
% use kernel simmetries
Kn_Gautschi = [Kn_Gautschi,...
               -conj(fliplr(Kn_Gautschi(1:qbar/4))),...
               -Kn_Gautschi(2:qbar/4+1),...
               conj(fliplr(Kn_Gautschi(2:qbar/4)))];
Kn_Gautschi = Kn_Gautschi.*dz;
Rq = zeros(1,p);
for ell = 1:p-1
  Rq_p = (1+sqrt(2))/qbar*max(abs(Kn_Gautschi*fellZA));
  Rq(ell:p) = Rq(ell:p)+Rq_p*normV(p:-1:ell)/2^(ell*s);
  Kn_Gautschi = Kn_Gautschi.*(z+1)/(2*ell);
end
Rq_p = (1+sqrt(2))/qbar*max(abs(Kn_Gautschi*fellZA));
Rq(p) = Rq(p)+Rq_p*normV(p)/2^(p*s);
if (any(Rq > tol))
  Rq = Inf(1,p);
end
