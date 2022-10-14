% Example of evolutionary advection--diffusion--reaction equation
% (see [CCZ22, Sec. 4.2])
%
% Equation:
% \partial_t u = \epsilon*\Delta u +
%                \alpha*(\partial_{x1} + \partial_{x2} + \partial_{x3})u +
%                g(u)
%
% u0 = 64*x1*(1-x1)*x2*(1-x2)*x3*(1-x3)
%
% where g is a nonlinear function of u.
% Spatial domain: [0,1]^3
% Time domain: [0,T] = [0,0.1]
% Boundary conditions: Homogeneous Dirichlet
% Space discretization: second order centered finite differences
% Time integration method: exponential Euler, ETD2RK
%
% [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%         A mu-mode approach for exponential integrators: actions of
%         phi-functions of Kronecker sums, Submitted, 2022

clear all
addpath('./integrators','./utils')
addpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
        './ext_routines/phipmsimuliom','./ext_routines/KronPACK')
fprintf('\n---- Evolutionary ADR ----\n')

routinerange = {'phiks','kiops','bamphi','phipmsimuliom'};
d = 3;
T = 0.1;
nrange = [20,64,81,100,121]; % temporal order on nrange(1)
                             % cpu time on nrange(2:end)
n = ones(1,d)*nrange(1);
epsilon = 0.5;
alpha = 10*ones(1,d);
[A,Kfun,g,gv,U0,U0v] = adr_function(epsilon,alpha,n);

method = 'expRK1s1_row';
order = 1;
tsrange = (300:100:700);

Uref = U0*exp(T);
fprintf('\nMethod: exponential Euler (order 1)\n')
counter_r = 0;
fprintf('\n## Check order of convergence ##')
for rou = routinerange
  routine = rou{1};
  fprintf(['\nRoutine: ',routine,'\n'])
  counter_r = counter_r + 1;
  counter = 0;
  for ts = tsrange
    fprintf('ts = %i\n',ts)
    counter = counter+1;
    if strcmp(routine,'phiks')
      U = feval([method,'_',routine],T,ts,A,U0,g);
    else
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
    end
      err(counter_r,counter) = norm(U(:)-Uref(:),inf)/norm(Uref(:),inf);
  end
end
sym = {'s','*','+','o'};
figure;
for i = 1:length(routinerange)
  loglog(tsrange,err(i,:),sym{i},'DisplayName',routinerange{i})
  hold on
end
loglog(tsrange,err(1,end)*(tsrange/tsrange(end)).^(-order),'--k',...
       'DisplayName',sprintf('order %i',order))
legend;
xlabel('time steps')
ylabel('error')
title('exponential Euler - Order of convergence')
drawnow

fprintf('\n## Wall-clock time ##\n')
ts = 250;
counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  counter_r = 0;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v,h] = adr_function(epsilon,alpha,n);
  for rou = routinerange
    routine = rou{1};
    fprintf(['n = %i, ',routine,'\n'],n(1))
    counter_r = counter_r+1;
    if strcmp(routine,'phiks')
      tic
      U = feval([method,'_',routine],T,ts,A,U0,g);
      cpu(counter_r,counter) = toc;
    else
      tic
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
      cpu(counter_r,counter) = toc;
    end
  end
end

figure;
for i = 1:length(routinerange)
  loglog(nrange(2:end),cpu(i,:),[sym{i},'-'],'DisplayName',routinerange{i})
  hold on
end
legend;
xlabel('n')
ylabel('wall-clock time')
title('exponential Euler - Wall-clock time')
drawnow
n = ones(1,d)*nrange(1);
[A,Kfun,g,gv,U0,U0v] = adr_function(epsilon,alpha,n);

method = 'expRK2s2_row';
order = 2;
tsrange = (200:50:400);

Uref = U0*exp(T);
fprintf('\nMethod: ETD2RK (order 2)\n')
counter_r = 0;
fprintf('\n## Check order of convergence ##')

for rou = routinerange
  routine = rou{1};
  fprintf(['\nRoutine: ',routine,'\n'])
  counter_r = counter_r + 1;
  counter = 0;
  for ts = tsrange
    fprintf('ts = %i\n',ts)
    counter = counter+1;
    if strcmp(routine,'phiks')
      U = feval([method,'_',routine],T,ts,A,U0,g);
    else
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
    end
    err(counter_r,counter) = norm(U(:)-Uref(:),inf)/norm(Uref(:),inf);
  end
end
sym = {'s','*','+','o'};
figure;
for i = 1:length(routinerange)
  loglog(tsrange,err(i,:),sym{i},'DisplayName',routinerange{i})
  hold on
end
loglog(tsrange,err(1,end)*(tsrange/tsrange(end)).^(-order),'--k',...
       'DisplayName',sprintf('order %i',order))
legend;
xlabel('time steps')
ylabel('error')
title('ETD2RK - Order of convergence')
drawnow

fprintf('\n## Wall-clock time ##\n')
ts = 100;
counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  counter_r = 0;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v] = adr_function(epsilon,alpha,n);
  for rou = routinerange
    routine = rou{1};
    fprintf(['n = %i, ',routine,'\n'],n(1))
    counter_r = counter_r+1;
    if strcmp(routine,'phiks')
      tic
      [U,scal,quadn,cost] = feval([method,'_',routine],T,ts,A,U0,g);
      cpu(counter_r,counter) = toc;
      scalLCP(counter) = mean(scal);
      quadnLCP(counter) = mean(quadn);
      costLCP(counter) = mean(cost);
    else
      tic
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
      cpu(counter_r,counter) = toc;
    end
  end
end

figure;
for i = 1:length(routinerange)
  loglog(nrange(2:end),cpu(i,:),[sym{i},'-'],'DisplayName',routinerange{i})
  hold on
end
legend;
xlabel('n')
ylabel('wall-clock time')
title('ETD2RK - Wall-clock time')
drawnow

fprintf('\n## Comparison phi-functions on the same vector ##\n')

counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v] = adr_function(epsilon,alpha,n);
  routine = 'phiks';
  method = 'expRK2s2_column';
  fprintf(['n = %i, ',routine,'\n'],n(1))
  counter_r = counter_r+1;
  tic
  [U,scal,quadn,cost] = feval([method,'_',routine],T,ts,A,U0,g);
  cpuPSV(counter) = toc;
  scalPSV(counter) = mean(scal);
  quadnPSV(counter) = mean(quadn);
  costPSV(counter) = mean(cost);
end
fprintf('\nETD2RK phiks (linear combination of phi-functions)\n')
fprintf('n = '),fprintf('%i ',nrange(2:end)),fprintf('\n');
fprintf('Wall-clock time = '),fprintf('%.2f ',cpu(1,:)),fprintf('\n');
fprintf('T = '),fprintf('%.1f ',costLCP),fprintf('\n');
fprintf('\nETD2RK phiks (phi-functions on the same vector)\n')
fprintf('n = '),fprintf('%i ',nrange(2:end)),fprintf('\n');
fprintf('Wall-clock time = '),fprintf('%.2f ',cpuPSV(1,:)),fprintf('\n');
fprintf('T = '),fprintf('%.1f ',costPSV),fprintf('\n');

rmpath('./integrators','./utils')
rmpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
       './ext_routines/phipmsimuliom','./ext_routines/KronPACK')
