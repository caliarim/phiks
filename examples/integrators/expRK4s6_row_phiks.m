function [U,s,q,cost] = expRK4s6_row_2_phiks(tstar,ts,A,U,g);
% Method: exponential Runge Kutta of order 4 (Luan)

  if ts == 20
    kappa = 2^(29); % brusselator
  else
    kappa = 2^0;
  end
  c(2)= 1/3;
  c(3)=1/3;
  c(4)=2*c(3);
  c(6)=1;
  c(5)=(4*c(6)-3)/(6*c(6)-4);
  tau = tstar / ts;
  tol = kappa*tau^5;
  t = 0;
  V{1}{1} = U{1};
  V{2}{1} = U{2};
  s = NaN(1,ts);
  q = NaN(1,ts);
  cost = NaN(1,ts);
  for j = 1:ts
    gtU_1 = g{1}(t,V{1}{1},V{2}{1});
    normU_1 = norm(V{1}{1}(:));
    V{1}{2} = c(2)*tau*gtU_1;
    [U2{1},s1_1,q1_1,c1_1] = phiks(c(2)*tau,A{1},V{1},1,tol*normU_1,1);
    gtU_2 = g{2}(t,V{1}{1},V{2}{1});
    normU_2 = norm(V{2}{1}(:));
    V{2}{2} = c(2)*tau*gtU_2;
    [U2{2},s1_2,q1_2,c1_2] = phiks(c(2)*tau,A{2},V{2},1,tol*normU_2,1);
    s1 = (s1_1+s1_2)/2;
    q1 = (q1_1+q1_2)/2;
    c1 = (c1_1+c1_2);
    D2_1 = g{1}(t+c(2)*tau,U2{1},U2{2})-gtU_1;
    V{1}{2} = c(4)*tau*gtU_1;
    V{1}{3} = c(4)^2*tau*D2_1/c(2);
    [U4{1},U3{1},s2_1,q2_1,c2_1] = phiks(c(4)*tau,A{1},V{1},2,tol*normU_1,2);
    D2_2 = g{2}(t+c(2)*tau,U2{1},U2{2})-gtU_2;
    V{2}{2} = c(4)*tau*gtU_2;
    V{2}{3} = c(4)^2*tau*D2_2/c(2);
    [U4{2},U3{2},s2_2,q2_2,c2_2] = phiks(c(4)*tau,A{2},V{2},2,tol*normU_2,2);
    s2 = (s2_1+s2_2)/2;
    q2 = (q2_1+q2_2)/2;
    c2 = (c2_1+c2_2);
    D3_1 = g{1}(t+c(3)*tau,U3{1},U3{2})-gtU_1;
    D4_1 = g{1}(t+c(4)*tau,U4{1},U4{2})-gtU_1;
    V{1}{2} = c(6)*tau*gtU_1;
    V{1}{3} = c(6)^2*tau*(-D3_1*c(4)/c(3)+D4_1*c(3)/c(4))/(c(3)-c(4));
    V{1}{4} = c(6)^3*tau*(D3_1/c(3)-D4_1/c(4))*2/(c(3)-c(4));
    [U6{1},U5{1},s3_1,q3_1,c3_1] = phiks(c(6)*tau,A{1},V{1},3,tol*normU_1,2);
    D3_2 = g{2}(t+c(3)*tau,U3{1},U3{2})-gtU_2;
    D4_2 = g{2}(t+c(4)*tau,U4{1},U4{2})-gtU_2;
    V{2}{2} = c(6)*tau*gtU_2;
    V{2}{3} = c(6)^2*tau*(-D3_2*c(4)/c(3)+D4_2*c(3)/c(4))/(c(3)-c(4));
    V{2}{4} = c(6)^3*tau*(D3_2/c(3)-D4_2/c(4))*2/(c(3)-c(4));
    [U6{2},U5{2},s3_2,q3_2,c3_2] = phiks(c(6)*tau,A{2},V{2},3,tol*normU_2,2);
    s3 = (s3_1+s3_2)/2;
    q3 = (q3_1+q3_2)/2;
    c3 = (c3_1+c3_2);
    D5_1 = g{1}(t+c(5)*tau,U5{1},U5{2})-gtU_1;
    D6_1 = g{1}(t+c(6)*tau,U6{1},U6{2})-gtU_1;
    V{1}{2} = tau*gtU_1;
    V{1}{3} = tau*(-D5_1*c(6)/c(5)+D6_1*c(5)/c(6))/(c(5)-c(6));
    V{1}{4} = tau*2*(D5_1/c(5)-D6_1/c(6))/(c(5)-c(6));
    [V{1}{1},s4_1,q4_1,c4_1] = phiks(tau,A{1},V{1},3,tol*normU_1,1);
    D5_2 = g{2}(t+c(5)*tau,U5{1},U5{2})-gtU_2;
    D6_2 = g{2}(t+c(6)*tau,U6{1},U6{2})-gtU_2;
    V{2}{2} = tau*gtU_2;
    V{2}{3} = tau*(-D5_2*c(6)/c(5)+D6_2*c(5)/c(6))/(c(5)-c(6));
    V{2}{4} = tau*2*(D5_2/c(5)-D6_2/c(6))/(c(5)-c(6));
    [V{2}{1},s4_2,q4_2,c4_2] = phiks(tau,A{2},V{2},3,tol*normU_2,1);
    s4 = (s4_1+s4_2)/2;
    q4 = (q4_1+q4_2)/2;
    c4 = (c4_1+c4_2);
    s(j) = (s1+s2+s3+s4)/4;
    q(j) = (q1+q2+q3+q4)/4;
    cost(j) = (c1+c2+c3+c4)/4;
    t = t+tau;
  end
  U{1} = V{1}{1};
  U{2} = V{2}{1};