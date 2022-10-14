% Example of Allen--Cahn equation
% (see [CCZ22, Sec. 4.3])
%
% Equation:
% \partial_t u = \Delta u +
%                (1/\epsilon^2)*u*(1-u^2)
%
% with "star shaped" initial datum u0 .
% Spatial domain: [0,1]^2
% Time domain: [0,T] = [0,0.025]
% Boundary conditions: Homogeneous Neumann
% Space discretization: second order centered finite differences
% Time integration method: third order exponential integrator
%
% [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%         A mu-mode approach for exponential integrators: actions of
%         phi-functions of Kronecker sums, Submitted, 2022

clear all
addpath('./integrators','./utils')
addpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
        './ext_routines/phipmsimuliom','./ext_routines/sylvphi',...
        './ext_routines/KronPACK')
fprintf('\n---- Allen--Cahn equation ----\n')

routinerange = {'phiks','kiops','bamphi','phipmsimuliom','sylvphi'};
d = 2;
T = 0.025;
nrange = [21,351,451,551,651]; % temporal order on nrange(1)
                               % cpu time on nrange(2:end)
n = ones(1,d)*nrange(1);
epsilon = 0.05;
N = 7;
alpha = 0.75;
[A,Kfun,g,gv,U0,U0v] = allencahn_function(epsilon,N,alpha,n,2);

method = 'expRK3s3_column';
order = 3;
tsrange = (100:25:200);
load('ref_allencahn_21.mat')

fprintf('\nMethod: exponential integrator of order 3 \n')
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
    if (strcmp(routine,'phiks') || strcmp(routine,'sylvphi'))
      U = feval([method,'_',routine],T,ts,A,U0,g);
    else
      U = feval([method,'_',routine],T,ts,Kfun,U0v,gv);
    end
    err(counter_r,counter) = norm(U(:)-Uref(:),inf)/norm(Uref(:),inf);
  end
end
sym = {'s','*','+','o','>'};
figure;
for i = 1:length(routinerange)
  loglog(tsrange,err(i,:),sym{i},'DisplayName',routinerange{i})
  hold on
end
loglog(tsrange,err(1,end)*(tsrange/tsrange(end)).^(-order),'--k','DisplayName',sprintf('order %i',order))
legend;
xlabel('time steps')
ylabel('error')
title('exponential integrator of order 3 - Order of convergence')
drawnow
fprintf('\n## Wall-clock time ##\n')
ts = 20;
counter = 0;
for n = nrange(2:end)
  counter = counter+1;
  counter_r = 0;
  n = ones(1,d)*n;
  [A,Kfun,g,gv,U0,U0v] = allencahn_function(epsilon,N,alpha,n,2);
  for rou = routinerange
    routine = rou{1};
    fprintf(['n = %i, ',routine,'\n'],n(1))
    counter_r = counter_r+1;
    if (strcmp(routine,'phiks') || strcmp(routine,'sylvphi'))
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
title('exponential integrator of order 3 - Wall-clock time')
drawnow

rmpath('./integrators','./utils')
rmpath('../','./ext_routines/bamphi','./ext_routines/kiops',...
        './ext_routines/phipmsimuliom','./ext_routines/sylvphi',...
        './ext_routines/KronPACK')
