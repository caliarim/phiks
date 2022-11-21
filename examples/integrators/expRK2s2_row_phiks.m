function [U,s,q,cost] = expRK2s2_row_phiks(tstar,ts,A,U,g)
% Method: ETD2RK

  kappa = 2^17; % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^3;
  V{1} = U;
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);
  for j = 1:ts
    gtU = g(t,V{1});
    V{2} = tau*gtU;
    normU = norm(V{1}(:));
    [U2,s1,q1,c1] = phiks(tau,A,V,1,tol*normU,1);
    V{3} = tau*(g(t+tau,U2)-gtU);
    [V{1},s2,q2,c2] = phiks(tau,A,V,2,tol*normU,1);
    s(j) = (s1+s2)/2;
    q(j) = (q1+q2)/2;
    cost(j) = c1+c2;
    t = t+tau;
  end
  U = V{1};
