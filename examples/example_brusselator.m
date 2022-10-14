% Example of Brusselator equation
% (see [CCZ22, Sec. 4.4])
%
% Equation:
% \partial_t u = d1*\Delta u - (b+1)*u + a + u^2*v
% \partial_t v = d2*\Delta v + b*u - u^2*v
% u0 = 64^2*x1^2*(1-x1)^2*x2^2*(1-x2)^2*x3^2*(1-x3)^2
% v0 = c
%
% Spatial domain: [0,1]^3
% Time domain: [0,T] = [0,1]
% Boundary conditions: Homogeneous Neumann
% Space discretization: fourth order centered finite differences
% Time integration method: fourth order exponential integrators
%
% [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%         A mu-mode approach for exponential integrators: actions of
%         phi-functions of Kronecker sums, Submitted, 2022

clear all
addpath('./integrators','./utils')
addpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
        './ext_routines/phipmsimuliom','./ext_routines/KronPACK')
fprintf('\n---- Brusselator equation ----\n')

routinerange = {'phiks','kiops','bamphi','phipmsimuliom'};
d = 3;
T = 1;
nrange = [11,40,60,80,100]; % temporal order on nrange(1)
                            % cpu time on nrange(2:end)
n = ones(1,d)*nrange(1);
d1 = 0.02;
d2 = 0.02;
c = 1;
a = 1;
b = 3;
[A,Kfun,g,gv,U0,U0v] = brusselator_function(a,b,c,d1,d2,n,4);
clear Kfun % create correct one inside loop

method = 'expRK4s5_mix';
order = 4;
tsrange = (40:10:80);
load('ref_brusselator_11_FD4.mat')

fprintf('\nMethod: Hochbruck--Ostermann of order 4 \n')
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
      temp = [U{1}(:);U{2}(:)];
    else
      tau = T/ts;
      KK{1} = tau*kronsum(A{1});
      KK{2} = tau*kronsum(A{2});
      prodn = prod(n);
      Kfun = @(u) [KK{1}*u(1:prodn);...
                   KK{2}*u(prodn+1:2*prodn)];
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
      temp = U(:);
    end
    tempref = [Uref{1}(:);Uref{2}(:)];
    err(counter_r,counter) = norm(temp(:)-tempref(:),inf)/norm(tempref(:),inf);
  end
end
sym = {'s','*','+','o'};
figure;
for i = 1:length(routinerange)
  loglog(tsrange,err(i,:),sym{i},'DisplayName',routinerange{i})
  hold on
end
loglog(tsrange,err(1,end)*(tsrange/tsrange(end)).^(-order),'--k','DisplayName',sprintf('order %i',order))
legend;
xlabel('time steps')
ylabel('error')
title('Hochbruck--Ostermann of order 4 - Order of convergence')
drawnow

fprintf('\n## Wall-clock time ##\n')
ts = 20;
counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  counter_r = 0;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v] = brusselator_function(a,b,c,d1,d2,n,4);
  clear Kfun % create correct one inside loop
  for rou = routinerange
    routine = rou{1};
    fprintf(['n = %i, ',routine,'\n'],n(1))
    counter_r = counter_r+1;
    if strcmp(routine,'phiks')
      tic
      U = feval([method,'_',routine],T,ts,A,U0,g);
      cpu(counter_r,counter) = toc;
    else
      tau = T/ts;
      KK{1} = tau*kronsum(A{1});
      KK{2} = tau*kronsum(A{2});
      prodn = prod(n);
      Kfun = @(u) [KK{1}*u(1:prodn);...
                   KK{2}*u(prodn+1:2*prodn)];
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
title('Hochbruck--Ostermann of order 4 - Wall-clock time')
drawnow

n = ones(1,d)*nrange(1);
[A,Kfun,g,gv,U0,U0v] = brusselator_function(a,b,c,d1,d2,n,4);

method = 'expRK4s6_row';
order = 4;
tsrange = (40:10:80);
load('ref_brusselator_11_FD4.mat')

fprintf('\nMethod: expRK4s6 (Luan) of order 4 \n')
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
      temp = [U{1}(:);U{2}(:)];
    else
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
      temp = U(:);
    end
    tempref = [Uref{1}(:);Uref{2}(:)];
    err(counter_r,counter) = norm(temp(:)-tempref(:),inf)/norm(tempref(:),inf);
  end
end
sym = {'s','*','+','o'};
figure;
for i = 1:length(routinerange)
  loglog(tsrange,err(i,:),sym{i},'DisplayName',routinerange{i})
  hold on
end
loglog(tsrange,err(1,end)*(tsrange/tsrange(end)).^(-order),'--k','DisplayName',sprintf('order %i',order))
legend;
xlabel('time steps')
ylabel('error')
title('expRK4s6 (Luan) of order 4 - Order of convergence')
drawnow

fprintf('\n## Wall-clock time ##\n')
ts = 20;
counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  counter_r = 0;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v] = brusselator_function(a,b,c,d1,d2,n,4);
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
title('expRK4s6 (Luan) of order 4 - Wall-clock time')
drawnow
rmpath('./integrators','./utils')
rmpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
        './ext_routines/phipmsimuliom','./ext_routines/KronPACK')
