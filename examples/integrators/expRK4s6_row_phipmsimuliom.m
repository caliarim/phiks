function [U,m_kry] = expRK4s6_row_phipmsimuliom(tstar,ts,Afun,U,g);
% Method: exponential Runge Kutta of order 4 (Luan)

  kappa = 2^(-8); %brusselator
  c(2)= 1/3;
  c(3)=1/3;
  c(4)=2*c(3);
  c(6)=1;
  c(5)=(4*c(6)-3)/(6*c(6)-4);
  tau = tstar / ts;
  tol = kappa*tau^5;
  t = 0;
  V(:,1) = U;
  m_kry1 = 1;
  m_kry2 = 1;
  m_kry3 = 1;
  m_kry4 = 1;
  iom = 2;
  for n = 1:ts
    gtU = g(t,V(:,1));
    normU = norm(V(:,1));
    krytol = tol*normU;
    V(:,2) = gtU;
    [U2,m_kry1] = phipm_simul_iom(c(2)*tau,Afun,V(:,1:2),krytol,m_kry1,iom);
    D2 = g(t+c(2)*tau,U2)-gtU;
    V(:,2) = gtU;
    V(:,3) = D2/c(2)/tau;
    [U34,m_kry2] = phipm_simul_iom([c(3),c(4)]*tau,Afun,V(:,1:3),krytol,...
                                   m_kry2,iom);
    D3 = g(t+c(3)*tau,U34(:,1))-gtU;
    D4 = g(t+c(4)*tau,U34(:,2))-gtU;
    V(:,2) = gtU;
    V(:,3) = ((-D3*c(4)/c(3)+D4*c(3)/c(4))/(c(3)-c(4)))/tau;
    V(:,4) = ((D3/c(3)-D4/c(4))*2/(c(3)-c(4)))/tau^2;
    [U56,m_kry3] = phipm_simul_iom([c(5),c(6)]*tau,Afun,V(:,1:4),...
                                   krytol,m_kry3,iom);
    D5 = g(t+c(5)*tau,U56(:,1))-gtU;
    D6 = g(t+c(6)*tau,U56(:,2))-gtU;
    V(:,2) = gtU;
    V(:,3) = ((-D5*c(6)/c(5)+D6*c(5)/c(6))/(c(5)-c(6)))/tau;
    V(:,4) = (2*(D5/c(5)-D6/c(6))/(c(5)-c(6)))/tau^2;
    [V(:,1),m_kry4] = phipm_simul_iom(tau,Afun,V(:,1:4),krytol,m_kry4,iom);
    t = t+tau;
  end
  U = V(:,1);
  m_kry = [m_kry1,m_kry2,m_kry3,m_kry4];
