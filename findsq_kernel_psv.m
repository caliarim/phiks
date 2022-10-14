function [sopt,qopt,costopt] = findsq_kernel_psv(GammaK,p,normV,tol,d,shat)
%
% function [sopt,qopt,costopt] = findsq_kernel_psv(GammaK,p,normV,tol,d,shat)
%
% Detect optimal scaling and number of quadrature points for
% phi-functions applied to the same vector.
% detect max scaling
rho = max(abs(GammaK));
smax = max(ceil(log2(rho)),0);
costopt = Inf;
qmax = 12;
qmin = 3;
r = 1.25;
qbar = 52; % 4*k
theta = linspace(0,2*pi,qbar+1);
theta = theta(1:qbar);
eitheta = exp(1i*theta);
c = 1; % ellipse with foci \pm c
a = r+c*c/(4*r);
b = r-c*c/(4*r);
z = r*eitheta+c*c/4/r./eitheta;
dz = (r*eitheta-c*c/4/r./eitheta); % 1i was removed by the absolute value
zbar = z(1:qbar/4+1);
[Z,A] = ndgrid(z,GammaK);
qtemp = qmax;
q = qmax;
s = shat-1;
smax = max(smax, s);
while ((s <= smax) && (q > qmin))
  Rq = Inf(1,p);
  fellZA = 1/2*exp((1-Z).*A/2/2^s);
  if (any(isinf(fellZA(:))))
    Rqtemp = Inf(1,p);
  else
    Rqtemp = esterr_kernel_psv(qtemp,p,GammaK,normV,tol,zbar,qbar,fellZA,z,dz);
  end
  while (all(~isinf(Rqtemp)) && qtemp > qmin)
    Rq = Rqtemp;
    q = qtemp;
    qtemp = qtemp - 1;
    Rqtemp = esterr_kernel_psv(qtemp,p,GammaK,normV,tol,zbar,qbar,fellZA,z,dz);
  end
  if (all(~isinf(Rqtemp)))
    Rq = Rqtemp;
    q = qtemp;
  end
  if (all(~isinf(Rq)))
    cost = q+p*s+shat;
    if (cost <= costopt)
      qopt = q;
      sopt = s;
      costopt = cost;
    else
      break
    end
  end
  s = s+1;
  GammaK = GammaK/2;
  tol = tol.*2.^(1:p);
end
if (isinf(costopt))
  sopt = smax;
  qopt = qmax;
  costopt = qopt+p*sopt+shat;
  warning(sprintf('Reached the maximum scaling and number of quadrature points.\n Try to reduce the matrix norm or to increase the tolerance.'))
end
