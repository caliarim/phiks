function [U,m_kry] = expRK3s3_column_phipmsimuliom(tstar,ts,Afun,U,g);
% Method: exponential Runge Kutta of order 3

  kappa = 2^(18); % allencahn
  c(2)=1/4;
  c(3)=1/2;
  gamma = (3*c(3)-2)*c(3)/(2-3*c(2))/c(2);
  tau = tstar / ts;
  tol = kappa*tau^4;
  t = 0;
  m_kry1 = 1;
  m_kry2 = 1;
  m_kry3 = 1;
  iom = 2;
  for j = 1:ts
    normU = norm(U(:));
    krytol = tol*normU;
    gtU = g(t,U);
    FtU(:,2) = Afun(U)+gtU;
    [aux1,m_kry1] = phipm_simul_iom(tau*[c(2),c(3),1],Afun,FtU,krytol,m_kry1,iom);
    U2 = U+aux1(:,1);
    D2(:,3) = g(t+c(2)*tau,U2)-gtU;
    [aux2,m_kry2] = phipm_simul_iom(tau*[c(2),c(3),1],Afun,D2/tau,krytol,m_kry2,iom);
    U3 = U+aux1(:,2)+(gamma/c(2))*aux2(:,1)+(1/c(2))*aux2(:,2);
    D3(:,3) = g(t+c(3)*tau,U3)-gtU;
    [aux3,m_kry3] = phipm_simul_iom(tau,Afun,D3/tau,krytol,m_kry3,iom);
    U = U+aux1(:,3)+gamma/(gamma*c(2)+c(3))*aux2(:,3)+1/(gamma*c(2)+c(3))*aux3;
    t = t+tau;
  end
  m_kry = [m_kry1,m_kry2,m_kry3];
