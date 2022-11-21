function [U,s,q,cost] = expRK2s2_column_phiks(tstar,ts,A,U,g)
%
% function U = expRK2s2_column(tstar,ts,A,U,g)
%
  %kappa = 2^10; % adr
  kappa = 2^13; % adr
  t = 0;
  tau = tstar/ts;
  tol = kappa*tau^3;
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);
  for j = 1:ts
    gtU = g(t,U);
    V = tau*(kronsumv(U,A)+gtU);
    normU = norm(U(:));
    [aux1,s1,q1,c1] = phiks(tau,A,V,1,tol*normU,1); % expEuler
    U2 = U+aux1{2};
    D2 = tau*(g(t+tau,U2)-gtU);
    [aux2,s2,q2,c2] = phiks(tau,A,D2,2,tol*normU,1);
    U = U2+aux2{3};
    s(j) = (s1+s2)/2;
    q(j) = (q1+q2)/2;
    cost(j) = c1+c2;
    t = t+tau;
  end
