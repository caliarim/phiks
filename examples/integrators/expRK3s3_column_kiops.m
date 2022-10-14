function [U,m_kry] = expRK3s3_column_kiops(tstar,ts,Afun,U,g);
% Method: exponential Runge Kutta of order 3

  kappa = 2^12; % allencahn
  c(2)=1/4;
  c(3)=1/2;
  gamma = (3*c(3)-2)*c(3)/(2-3*c(2))/c(2);
  tau = tstar / ts;
  tol = kappa*tau^4;
  t = 0;
  m_kry1 = [];
  m_kry2 = [];
  m_kry3 = [];
  for j = 1:ts
    normU = norm(U(:));
    krytol = tol*normU;
    gtU = g(t,U);
    FtU(:,2) = Afun(U)+gtU;
    [aux1,m_kry1] = kiops(tau*[c(2),c(3),1],Afun,FtU,krytol,m_kry1);
    U2 = U+c(2)*tau*aux1(:,1);
    D2(:,3) = g(t+c(2)*tau,U2)-gtU;
    [aux2,m_kry2] = kiops(tau*[c(2),c(3),1],Afun,D2,krytol,m_kry2);
    U3 = U+tau*(c(3)*aux1(:,2)+c(2)*gamma*aux2(:,1)+...
                     c(3)^2/c(2)*aux2(:,2));
    D3(:,3) = g(t+c(3)*tau,U3)-gtU;
    [aux3,m_kry3] = kiops(tau,Afun,D3,krytol,m_kry3);
    U = U+tau*(aux1(:,3)+gamma/(gamma*c(2)+c(3))*aux2(:,3)+...
               1/(gamma*c(2)+c(3))*aux3);
    t = t+tau;
  end
  m_kry = [m_kry1,m_kry2,m_kry3];
