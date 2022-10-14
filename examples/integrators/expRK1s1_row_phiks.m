function [U,s,q,cost] = expRK1s1_row_phiks(tstar,ts,A,U,g)
% Method: exponential Euler

  kappa = 2^11; % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^2;
  V{1} = U;
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);
  for j = 1:ts
    gtU = g(t,V{1});
    V{2} = tau*gtU;
    normU = norm(V{1}(:));
    [V{1},s(j),q(j),cost(j)] = phiks(tau,A,V,1,tol*normU);
    t = t+tau;
  end
  U = V{1};
