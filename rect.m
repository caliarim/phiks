function extreigs = rect (A, method)
%
% function extreigs = rect (A, method)
%
% A is assumed to be full. method == 'norm2' is norm2 (default),
% method == 'gersh' is Gershgorin disks. On output,  extreigs = [SR,LR,SI,LI].
if (nargin == 1)
  method = 'norm2';
end
n = length(A);
switch method
  case {'norm2'}
    mu = trace(A)/n;
    A(1:n+1:n*n) = A(1:n+1:n*n)-mu;
    SRLR = norm(A+A')/2;
    extreigs(1) = -SRLR+real(mu);
    extreigs(2) = SRLR+real(mu);
    SILI = norm(A-A')/2;
    extreigs(3) = -SILI+imag(mu);
    extreigs(4) = SILI+imag(mu);
  case {'gersh'}
    radiusR = zeros(n,1);
    radiusI = zeros(n,1);
    for i = 1:n
      radiusR(i) = sum(abs(A(i,:)'+A(:,i)))/2;
      radiusI(i) = sum(abs(A(i,:)'-A(:,i)))/2;
    end
    d = diag(A);
    rd = real(d);
    id = imag(d);
    extreigs = [min(rd - (radiusR - abs(rd))),...
                max(rd + (radiusR - abs(rd))),...
                min(id - (radiusI - abs(id))),...
                max(id + (radiusI - abs(id)))];
  otherwise
    error('rect: method not defined')
end
%!demo
%! A = randn(10)+1i*randn(10);
%! lambda = eig(A);
%! plot(real(lambda),imag(lambda),'*')
%! extreigs = rect(A,'gersh');
%! hold on
%! plot(extreigs([1,2,2,1,1]),extreigs([3,3,4,4,3]))
%! extreigs = rect(A,'norm2');
%! plot(extreigs([1,2,2,1,1]),extreigs([3,3,4,4,3]))
%! legend('eigenvalues','Gershgorin','norm2')
%!demo
%! A = toeplitz([-2,1,zeros(1,8)]);
%! lambda = eig(A);
%! plot(real(lambda),imag(lambda),'*')
%! extreigs = rect(A,'gersh');
%! hold on
%! plot(extreigs([1,2,2,1,1]),extreigs([3,3,4,4,3]))
%! extreigs = rect(A,'norm2');
%! plot(extreigs([1,2,2,1,1]),extreigs([3,3,4,4,3]))
%! legend('eigenvalues','Gershgorin','norm2')
