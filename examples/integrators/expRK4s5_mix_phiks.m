function [U,s,q,cost] = expRK4s5_mix_phiks(tstar,ts,A,U,g);
% Method: exponential Runge Kutta of order 4

  if ts == 20
    kappa = 2^22; % brusselator
  else
    kappa = 1;
  end
  c(2)=1/2;
  c(3)=1/2;
  c(4)=1;
  c(5)=1/2;
  tau = tstar / ts;
  tol = kappa*tau^5;
  t = 0;
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);

  for j = 1:ts
    normU_1 = norm(U{1}(:));
    normU_2 = norm(U{2}(:));

    gtU_1 = g{1}(t,U{1},U{2});
    gtU_2 = g{2}(t,U{1},U{2});

    FtU_1 = kronsumv(U{1},A{1})+gtU_1;
    FtU_2 = kronsumv(U{2},A{2})+gtU_2;

    [aux1_1,aux1half_1,s1_1,q1_1,c1_1] = phiks(tau,A{1},FtU_1,1,tol*normU_1,2);
    [aux1_2,aux1half_2,s1_2,q1_2,c1_2] = phiks(tau,A{2},FtU_2,1,tol*normU_2,2);

    s1 = (s1_1+s1_2)/2;
    q1 = (q1_1+q1_2)/2;
    c1 = c1_1+c1_2;

    U2_1 = U{1}+c(2)*tau*aux1half_1{2};
    U2_2 = U{2}+c(2)*tau*aux1half_2{2};
    D2_1 = g{1}(t+c(2)*tau,U2_1,U2_2)-gtU_1;
    D2_2 = g{2}(t+c(2)*tau,U2_1,U2_2)-gtU_2;

    [aux2_1,aux2half_1,s2_1,q2_1,c2_1] = phiks(tau,A{1},D2_1,2,tol*normU_1,2);
    [aux2_2,aux2half_2,s2_2,q2_2,c2_2] = phiks(tau,A{2},D2_2,2,tol*normU_2,2);
    s2 = (s2_1+s2_2)/2;
    q2 = (q2_1+q2_2)/2;
    c2 = c2_1+c2_2;

    U3_1 = U{1}+c(3)*tau*aux1half_1{2}+tau*aux2half_1{3};
    U3_2 = U{2}+c(3)*tau*aux1half_2{2}+tau*aux2half_2{3};
    D3_1 = g{1}(t+c(3)*tau,U3_1,U3_2)-gtU_1;
    D3_2 = g{2}(t+c(3)*tau,U3_1,U3_2)-gtU_2;

    [aux3_1,aux3half_1,s3_1,q3_1,c3_1] = phiks(tau,A{1},D3_1,2,tol*normU_1,2);
    [aux3_2,aux3half_2,s3_2,q3_2,c3_2] = phiks(tau,A{2},D3_2,2,tol*normU_2,2);
    s3 = (s3_1+s3_2)/2;
    q3 = (q3_1+q3_2)/2;
    c3 = c3_1+c3_2;

    U4_1 = U{1}+c(4)*tau*aux1_1{2}+tau*aux2_1{3}+tau*aux3_1{3};
    U4_2 = U{2}+c(4)*tau*aux1_2{2}+tau*aux2_2{3}+tau*aux3_2{3};
    D4_1 = g{1}(t+c(4)*tau,U4_1,U4_2)-gtU_1;
    D4_2 = g{2}(t+c(4)*tau,U4_1,U4_2)-gtU_2;

    [aux4_1,aux4half_1,s4_1,q4_1,c4_1] = phiks(tau,A{1},...
                                         {0,0,-1/4*D4_1,-D2_1-D3_1+D4_1},...
                                         3,tol*normU_1,2);
    [aux4_2,aux4half_2,s4_2,q4_2,c4_2] = phiks(tau,A{2},...
                                         {0,0,-1/4*D4_2,-D2_2-D3_2+D4_2},...
                                         3,tol*normU_2,2);
    s4 = (s4_1+s4_2)/2;
    q4 = (q4_1+q4_2)/2;
    c4 = c4_1+c4_2;

    U5_1 = U{1}+c(5)*tau*aux1half_1{2}+...
         tau*(1/2*aux2half_1{3}+1/2*aux3half_1{3}+4*aux4half_1+...
            1/4*aux2_1{3}+1/4*aux3_1{3}+aux4_1);
    U5_2 = U{2}+c(5)*tau*aux1half_2{2}+...
         tau*(1/2*aux2half_2{3}+1/2*aux3half_2{3}+4*aux4half_2+...
            1/4*aux2_2{3}+1/4*aux3_2{3}+aux4_2);
    D5_1 = g{1}(t+c(5)*tau,U5_1,U5_2)-gtU_1;
    D5_2 = g{2}(t+c(5)*tau,U5_1,U5_2)-gtU_2;

    [aux5_1,s5_1,q5_1,c5_1] = phiks(tau,A{1},...
                                    {0,0,-D4_1+4*D5_1,4*D4_1-8*D5_1},...
                                    3,tol*normU_1,1);
    [aux5_2,s5_2,q5_2,c5_2] = phiks(tau,A{2},...
                                    {0,0,-D4_2+4*D5_2,4*D4_2-8*D5_2},...
                                    3,tol*normU_2,1);
    s5 = (s5_1+s5_2)/2;
    q5 = (q5_1+q5_2)/2;
    c5 = c5_1+c5_2;

    U{1} = U{1}+tau*(aux1_1{2}+aux5_1);
    U{2} = U{2}+tau*(aux1_2{2}+aux5_2);
    s(j) = (s1+s2+s3+s4+s5)/5;
    q(j) = (q1+q2+q3+q4+q5)/5;
    cost(j) = c1+c2+c3+c4+c5;
    t = t+tau;
  end
