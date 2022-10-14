% Code validation for PHIKS (see [CCZ22, Sec. 4.1])
%
% [CCZ22] M. Caliari, F. Cassini, and F. Zivcovich,
%         A mu-mode approach for exponential integrators: actions of
%         phi-functions of Kronecker sums, Submitted, 2022

clear all
addpath('../','./ext_routines/kiops','./ext_routines/KronPACK')
fprintf('---- Code validation ----\n')

tol = 2^(-53);
p = 5;
tau = 1;
drange = [3,6];
nrange = [64,81,100,121;8,9,10,11];
counter_d = 0;

for d = drange
  counter_d = counter_d + 1;
  counter = 0;
  for n = nrange(counter_d,:)
    fprintf('\nd = %i, n = %i\n',d,n)
    counter = counter+1;
    for mu = 1:d
      x{mu} = linspace(0,1,n+2);
      h(mu) = 1/(n+1);
      x{mu} = x{mu}(2:n+1);
      A{mu} = (1+1i)/100*spdiags(ones(n,1)*[1,-2,1]/h(mu)^2,-1:1,n,n);
      fA{mu} = full(A{mu});
    end
    [X{1:d}] = ndgrid(x{1:d});
    V0 = 4096*(1+1i)*X{1}.*(1-X{1});
    for mu = 2:d
      V0 = V0.*X{mu}.*(1-X{mu});
    end
    [res1,res12,s,q,c] = phiks(tau,fA,V0,p,tol,2,1);
    K = kronsum(A);
    for ell = 0:p
      fprintf('ell = %i\n',ell)
      V = zeros(n^d,ell+1);
      V(:,ell+1) = V0(:);
      ref = kiops([tau/2,tau],K,V(:,1:ell+1),tol);
      err_psi(ell+1,counter,counter_d) = norm(res1{ell+1}(:)-ref(:,2),inf)/norm(ref(:,2),inf);
      err_psi12(ell+1,counter,counter_d) = norm(res12{ell+1}(:)-ref(:,1),inf)/norm(ref(:,1),inf);
    end
    srange_PSV(counter,counter_d) = s;
    qrange_PSV(counter,counter_d) = q;
    crange_PSV(counter,counter_d) = c;
    clear V
    V{1} = 0;
    for ell = 1:p
      V{ell+1} = V0;
    end
    [res2,res12,s,q,c] = phiks(tau,fA,V,p,tol,2,1);
    clear V;
    V(:,1) = zeros(n^d,1);
    for ell = 1:p
      V(:,ell+1) = V0(:)/tau^(ell);
    end
    ref = kiops([tau/2,tau],K,V(:,1:p+1),tol,[],[],[],false);
    err_psi_lc(counter,counter_d) = norm(res2(:)-ref(:,2),inf)/norm(ref(:,2),inf);
    err_psi_lc12(counter,counter_d) = norm(res12(:)-ref(:,1),inf)/norm(ref(:,1),inf);
    srange_LCP(counter,counter_d) = s;
    qrange_LCP(counter,counter_d) = q;
    crange_LCP(counter,counter_d) = c;
  end
  fprintf('\nd = %i: s, q and number of Tucker operators (phi-functions on the same vector)\n',d)
  fprintf('s = '),fprintf('%i ',srange_PSV(:,counter_d)),fprintf('\n');
  fprintf('q = '),fprintf('%i ',qrange_PSV(:,counter_d)),fprintf('\n');
  fprintf('T = '),fprintf('%i ',crange_PSV(:,counter_d)),fprintf('\n');
  fprintf('\nd = %i: s, q and number of Tucker operators (linear combination of phi-functions)\n',d)
  fprintf('s = '),fprintf('%i ',srange_LCP(:,counter_d)),fprintf('\n');
  fprintf('q = '),fprintf('%i ',qrange_LCP(:,counter_d)),fprintf('\n');
  fprintf('T = '),fprintf('%i ',crange_LCP(:,counter_d)),fprintf('\n');
  figure;
  for ell = 2:p+1
    semilogy(nrange(counter_d,:),err_psi(ell,:,counter_d),'*-','DisplayName',['ell = ',num2str(ell-1)])
    hold on
  end
  semilogy(nrange(counter_d,:),err_psi_lc(:,counter_d),'*-','DisplayName','linear combination')
  legend;
  xlabel('n')
  ylabel('difference')
  ylim([2e-16,2e-11])
  title(sprintf('d = %d, p = %d',d,p))
  drawnow
  figure;
  for ell = 2:p+1
    semilogy(nrange(counter_d,:),err_psi12(ell,:,counter_d),'*-','DisplayName',['ell = ',num2str(ell-1)])
    hold on
  end
  semilogy(nrange(counter_d,:),err_psi_lc12(:,counter_d),'*-','DisplayName','linear combination')
  legend;
  xlabel('n')
  ylabel('difference')
  ylim([2e-16,2e-11])
  title(sprintf('d = %d, p = %d, half scale',d,p))
  drawnow
end
rmpath('../','./ext_routines/kiops','./ext_routines/KronPACK')
