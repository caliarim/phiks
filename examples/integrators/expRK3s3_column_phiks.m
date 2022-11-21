function [U,s,q,cost] = expRK3s3_column_phiks(tstar,ts,A,U,g);
% Method: exponential Runge Kutta of order 3

  kappa = 2^(36); % allencahn
  c(2)=1/4;
  c(3)=1/2;
  gamma = (3*c(3)-2)*c(3)/(2-3*c(2))/c(2);
  tau = tstar / ts;
  tol = kappa*tau^4;
  t = 0;
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);
  for j = 1:ts
    normU = norm(U(:));
    gtU = g(t,U);
    FtU = kronsumv(U,A)+gtU;
    [aux1,aux1half,aux1quarter,s1,q1,c1] = phiks(tau,A,FtU,1,tol*normU,3);
    U2 = U+c(2)*tau*aux1quarter{2};
    D2 = g(t+c(2)*tau,U2)-gtU;
    [aux2,aux2half,aux2quarter,s2,q2,c2] = phiks(tau,A,D2,2,tol*normU,3);
    U3 = U+tau*(c(3)*aux1half{2}+c(2)*gamma*aux2quarter{3}+...
                     c(3)^2/c(2)*aux2half{3});
    D3 = g(t+c(3)*tau,U3)-gtU;
    [aux3,s3,q3,c3] = phiks(tau,A,D3,2,tol*normU,1);
    U = U+tau*(aux1{2}+gamma/(gamma*c(2)+c(3))*aux2{3}+...
               1/(gamma*c(2)+c(3))*aux3{3});
    s(j) = (s1+s2+s3)/3;
    q(j) = (q1+q2+q3)/3;
    cost(j) = c1+c2+c3;
    t = t+tau;
  end
