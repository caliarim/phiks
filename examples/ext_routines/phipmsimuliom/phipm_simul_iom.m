function [wout, m, stats] = phipm_simul_iom(t, Afun, u, tol, m, iom)
% Usage: [wout, m, stats] = phipm_simul_iom(t, Afun, u, tol, m, iom)
%
% Evaluates a linear combination of the phi functions evaluated at
% t*A acting on a set of input vectors.
%
% The evaluation expresses eveything in terms of the highest order
% phi function and evaluates the action of this on a vector using a
% Krylov technique and then computes w using the recurrence
% relation. 
%
% The size of the Krylov subspace is changed dynamically during the
% integration. The Krylov subspace is computed using the Arnoldi process.
%
% This file is a modified version of phipm_iom.m, developed by
% S. Gadreault and J.A. Pudykiewicz (JCP 2016).
% Similarly, phipm_iom.m is a modified version of phipm.m,
% developed by J. Niesen and W.M. Wright (ACM TOMS 2012).
%
% Modifications performed by Daniel R. Reynolds and Vu Luan, 
% Mathematics Department, Southern Methodist University 
% Fall 2017
%
% INPUTS:
%   t    - array of time values.  We assume that these are distinct
%          and in order, i.e. 0 < |t(1)| < |t(2)| < ... < |t(end)|
%   Afun - the action of the matrix argument of the phi functions
%   u    - the matrix with columns representing the vectors to be
%          multiplied by the phi functions
%   tol  - the convergence tolarance required
%   m    - an estimate of the appropriate Krylov size
%   iom  - orthogonalization length
%
% OUTPUTS:
%   wout - matrix of outputs corresponding to the desired linear
%           combination evaluated at the set of times specified
%           by the input array t:
%               w(:,1) = \sum_{l=0}^{p} t(1)^l * \phi_l(t(1)*A)*u(:,l) 
%               w(:,2) = \sum_{l=0}^{p} t(2)^l * \phi_l(t(2)*A)*u(:,l) 
%               ...
%               w(:,N) = \sum_{l=0}^{p} t(N)^l * \phi_l(t(N)*A)*u(:,l) 
%           where p is the number of columns from the input u, N is
%           the number of times supplied from the input t, and
%           \phi_l() is the l-th phi function. 
%   m     - size of the Krylov subspace corresponding to the
%           solution (in the final substep only).
%   stats - [# substeps, # rejected steps, # Krylov steps, # matrix exponentials]
%
% 14th October 2022
% This file has been slightly modified to accept in input a function that
% performs the action of the matrix [CCZ22].

% retain reusable internal structures in memory between calls to this function
persistent V int nnze

% n is the size of the original problem
% p is the number of vectors to multiply (and the number of phi functions)
%   (if we need to compute phi_0, phi_1 and phi_2 then p = 3)
[n, p] = size(u);

% nt is the number of time outputs
N = length(t);

% ensure that the first entry in t is nonzero
if (abs(t(1)) == 0)
   error('phipm_simul_iom: error, all entries of t must be nonzero');
end


% initialize output 
wout = zeros(n,N);

% Add extra column of zeros if p=1
if p == 1
   p = 2;
   u = [u, zeros(size(u))];
end

% Krylov parameters
mmax = 100;
mnew = m;

% Preallocate matrices V and int
if isempty(V) && isempty(int)
   V = zeros(n, mmax+1);
   int = zeros(n, p);
elseif numel(V) ~= n*(mmax+1)
   V = zeros(n, mmax+1);
   int = zeros(n, p);
elseif numel(int) ~= n*p
   int = zeros(n, p);
end

% determine number of nonzeros in A (for cost estimate)
if (isnumeric(Afun))
  if isempty(nnze)
    nnze = nnz(Afun);
  end
  Afun = @(x) Afun*x;
else
  nnze = 7*n;
end

% Initializing the variables for output statistics 
step = 0;
reject = 0;
krystep = 0;
exps = 0;

% Setting the safety factors and tolerance requirements
gamma = 0.8;
delta = 1.4;

% Temporary vectors used for toeplitz trick
cidx = (0:p-1)';
ridx = p:-1:1;
idx = cidx(:,ones(p,1)) + ridx(ones(p,1),:);

% initialize additional temporary variables
ireject = 0;         % flag indicating number of rejected steps
happy = false;       % flag indicating a 'happy' breakdown 
sgn = sign(t(end));  % time 'direction' (positive or negative)
tnow = 0;            % internal time [0,tout]
j = 0;               % counter for Krylov vector

% Set initial condition
w = u(:, 1);

% Initialize extra constants used for adaptivity
oldm = NaN; 
oldtau = NaN; 
omega = NaN;
orderold = true; 
kestold = true;

% Set initial timestep
tau = abs(t(1));

% loop through array of output times
for itout = 1:N

   % set current output time
   tout = abs(t(itout));

   % Iterate until we reach the desired output time 
   while (tnow < tout)

      % do not overstep tout
      tau = min(tout-tnow, tau);
      
      % Compute necessary starting information
      if (j == 0)

         % Initialize the matrices H, x and up
         % Note:  x(idx) is a lower triangular matrix with 1 on the
         % diagonal and (tnow^j / j!) on the j-th subdiagonal, so
         %    up(:,k) = sum_{j=0}^{p-k} (tnow^j / j! * u(:,k+j))
         H = zeros(mmax+p, mmax+p);
         x = [ zeros(1,p-1), cumprod([1,tnow./(1:p-1)]) ];
         up = u*x(idx);  % even if u has zero columns, the 
                         % subdiagonal can fill those colums in up

         % Compute the update factors
         % (the w_j vectors in Section 3.3 of the paper)
         int(:, 1) = w;
         for i = 1:p-1
            if (norm(int(:,i)) > 0)
               int(:,i+1) = Afun(int(:,i)) + up(:,i+1);
            else
               int(:,i+1) = up(:,i+1);
            end
         end

         % Compute normalization factor
         beta = norm(int(:,end));
         if (beta == 0)
            
            % Multiplying with a zero vector, hence result is zero
            % Finish all in one step
            reject = reject + ireject;
            step = step + 1;
            tau = tout - tnow;
            w = w + int(:,2:p-1)*cumprod(tau*1./(1:p-2)');
            break;

         end

         % The first Krylov basis vector
         V(:,1) = int(:,end)./beta;
         
      end

      % begin Arnoldi process
      while (j < m)

         % update column index, and compute matrix-vector product
         j = j+1;
         if (norm(V(:,j)) > 0)
            vv = Afun(V(:,j));
         else
            vv = V(:,j);
         end

         % update H and vv via incomplete orthogonalization;
         % increment krystep and recalculate s
         for i = max(1,j-iom):j
            H(i,j) = V(:,i)'*vv;
            vv = vv - H(i,j)*V(:,i);
         end
         krystep = krystep + 1;
         s = norm(vv);
         
         % Happy breakdown
         if (s < tol)
            happy = true;
            tau = tout - tnow;
            break;
         end
         
         % update H and V for next iteration
         H(j+1,j) = s;
         V(:,j+1) = vv./s;
         
      end

      % Keep a record of H, in case of failure
      H2 = H;

      % We use the vector e1 in the computations
      H(1,j+1) = 1;

      % Construct the augmented H matrix, store h for later
      for i = 1:p-1
         H(j+i,j+i+1) = 1;
      end
      h = H(j+1,j);
      H(j+1,j) = 0;
      
      % Compute the exponential of the augmented matrix, and increment counter
      [F,hnorm] = expmnorm( sgn*tau*H(1:j+p,1:j+p) );
      exps = exps+1;

      % Estimate local truncation error for internal step
      err = abs(beta*h*F(j,j+p));

      % Compute error per unit step
      oldomega = omega;
      omega = tout*err/(tau*tol);
              
      % Estimate temporal order of accuracy
      if ( (m == oldm) && (tau ~= oldtau) && (ireject >= 1) )
         order = max(1, log(omega/oldomega)/log(tau/oldtau));
         orderold = false;
      elseif ( orderold || (ireject == 0) )
         orderold = true;
         order = j/4;
      else
         orderold = true;
      end
      
      % Estimate convergence factor, k, for varying Krylov subspace size
      if ( (m ~= oldm) && (tau == oldtau) && (ireject >= 1) )
         kest = max(1.1, (omega/oldomega)^(1/(oldm-m)));
         kestold = false;
      elseif ( kestold || (ireject == 0) )
         kestold = true;
         kest = 2;
      else
         kestold = true;
      end

      % This if block is the main difference between fixed and variable m
      oldtau = tau; 
      oldm = m;
      if happy

         % Happy breakdown; wrap up
         omega = 0;
         taunew = tau;
         mnew = m;

      elseif ( (j == mmax) && (omega > delta) )
         
         % Krylov subspace too small and stepsize too large
         taunew = tau*(omega/gamma)^(-1/order);
         
      else

         % Determine optimal tau and m
         tauopt = tau*(omega/gamma)^(-1/order);
         mopt = max(1, ceil(j+log(omega/gamma)/log(kest)));
         nom = 5 + max(log(hnorm), 0)/log(2);  % number of mult's in expm
         
         % Cost [flops] of Arnoldi
         %cost = @(M,T) ((M+p)*nnze+(M^2+3*p+2)*n+nom*(M+p-1)^3)*ceil((tout-tnow)/T);

         % Cost [flops] of incomplete orthognalization is almost like Lanczos ...
         cost = @(M,T) ((M+p)*nnze+3*(M+p)*n+nom*(M+p-1)^3)*ceil((tout-tnow)/T);

         % Determine whether to vary tau or m
         if (cost(j, tauopt) < cost(mopt, tau))
            taunew = tauopt;
            mnew = m;
         else
            mnew = mopt;
            taunew = tau;
         end
         
      end
      
      % Check error against target
      if (omega <= delta)

         % Achieved the required tolerance; update solution:
         %    up = w+sum_{k=1}^{p-2} tau^k / k! * int(:,k+1)
         reject = reject + ireject;
         step = step + 1;
         up = w + int(:,2:p-1)*cumprod(tau*1./(1:p-2)');
         
         % Update solution using the corrected quantity
         F(j+1,j+p-1) = h*F(j,j+p);
         w = up + beta*V(:,1:j+1)*F(1:j+1,j+p-1);
         
         % Update t and counters
         tnow = tnow + tau;
         j = 0;
         ireject = 0;
         
      else

         % Did not achieve target tolerance, try again
         H = H2;
         ireject = ireject + 1;
         
      end
      
      % Safety factors for tau
      tau = max(tau/5, min(2*tau, taunew));
      
      % Safety factors for m
      m = max(1, min(mmax, max(floor(3/4*m), min(mnew, ceil(4/3*m)))));
      
   end

   % fill output with current solution
   wout(:,itout) = w;
   
end

% fill output statistics array
stats = [step, reject, krystep, exps];

% end of function
