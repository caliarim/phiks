function [ f, info ] = bamphi( t, A, At, u, opts, info )
%
% BAMPHI: Backward-stable Action of Matrix PHI-functions
% Version of Jan the 6th, 2022
%
% BAMPHI produces two outputs:
%
%   * f, which is the matrix s.t.
%
%        f( :,i ) = \sum_{j=0}^{p-1} t(i)^j * \phi_j(t*A) * u( :,j+1 )     ( 1 )
%
%   * info, that is a structure (not meant to be read by the user) with 4 fields:
%        -        A, information about the matrix A and the approximation params
%        - scal_ref, information about the timesteps-size BAMPHI worked with
%
% BAMPHI accepts the following inputs:
%
%   *    t: (mandatory) vector of timesteps as in ( 1 )
%
%   *    A: (mandatory) function handle such that  A( x ) = A  * x, with A as in ( 1 )
%
%   *   At: (optional/not recommended) function handle such that At( x ) = A' * x
%
%   *    u: (mandatory) n x p matrix as in ( 1 )
%
%   * opts: (optional)  is a structure whose fields are the desired options (see below)
%
%   * info: (optional)  info contains important information and it contains the
%                       suggested interpolation set. Suggestion: if you want to
%                       mess things up do it right and run bamphi once to learn
%                       how info is structured.
%
% BAMPHI accepts the following options:
%
%   * Error options
%        - opts.tol:      Tolerance set by the user. Default is 2^-53.
%        - opts.error:    Absolute or Relative error. Default is 'abs', other
%                         choice: 'rel'.
%        - opts.norm:     How the error is measured when checking for early
%                         termination. Default is inf, other choices are those
%                         accepted by Matlab's norm function.
%
%   * Technical options
%        - opts.low_stor: Set it to true if the problem is large enough for
%                         memory latencies to appear when storing several vectors
%                         in forming f. Default is false.
%        - opts.iom:      Set it to true for Incomplete Ortogonalization Method
%                         in Arnoldi process. Default is true.
%        - opts.scal_ref: Set it to true to allow scaling refinement. Default is
%                         true.
%
%   * Approximation options
%        - opts.DEG_MAX:  Maximum approximation polynomial's degree. Default is 128.
%        - opts.m_r_arn:  Maximum size of Krylov space. Default is min( 64, opts.DEG_MAX ).
%        - opts.r_arn:    Suggested size of Krylov space. Default is min( opts.m_r_arn, 64 ).
%        - opts.r_lan:    Maximum size of Krylov space when using Lanczos. Default is 64.
%
%   * Field of Values (FoV) options
%        - opts.skew:     Set to 1 if A is hermitian, -1 if skew Hermitian, 0
%                         else or unknown. Default is 0.
%        - opts.fov_tol:  Tolerance in computing FoV. Default is 1e-3.
%

  %% 01. CHECK INPUTS
  if ( nargin < 4 ), error('Not enough input arguments.'); end
  if ( nargin < 5 ) || ~isstruct( opts ), opts = []; end
  if ( nargin < 6 ) || ~isstruct( info ), info = []; end
  if all( t == 0 ),
    f = repmat( u( :,1 ), 1, length( t ) );
    return
  end
  [ n, p ] = size( u ); p = p - 1;

  if ( p > 0 )
    n_u = max( 1, 2^ceil( log2( sqrt( norm( u' * u ) ) ) ) ); % u'*u is a p x p matrix so it's faster than norm( u )
  else
    n_u = 1;
  end
  if ( n_u == 0 ),
    f = repmat( u( :,1 ), 1, length( t ) );
    return
  elseif ( n_u > 1 )
    u = u / n_u;
  end

  %% 02. CHECK OPTIONS
  % error options
  if ~isfield( opts,'tol'      ) || isempty( opts.tol      ), opts.tol      = 2^-53;                   end
  if ~isfield( opts,'error'    ) || isempty( opts.error    ), opts.error    = true;                    else, opts.error = strcmp( opts.error, 'rel' ); end
  if ~isfield( opts,'norm'     ) || isempty( opts.norm     ), opts.norm     = inf;                     end
  % technical options
  if ~isfield( opts,'low_stor' ) || isempty( opts.low_stor ), opts.low_stor = false;                   end
  if ~isfield( opts,'iom'      ) || isempty( opts.iom      ), opts.iom      = true;                    end
  if ~isfield( opts,'earl_ter' ) || isempty( opts.earl_ter ), opts.earl_ter = true;                    end
  if ~isfield( opts,'scal_ref' ) || isempty( opts.scal_ref ), opts.scal_ref = opts.earl_ter;           end
  if ~isfield( opts,'SCRF_LOW' ) || isempty( opts.SCRF_LOW ),
    if ~isfield( info,'scal_ref' ) || isempty( info.scal_ref ) || ~isfield( info.scal_ref,'SCRF_LOW' ) || isempty( info.scal_ref.SCRF_LOW ),
      opts.SCRF_LOW = 0.60;
    else
      opts.SCRF_LOW = info.scal_ref.SCRF_LOW;
    end
  end
  if ~isfield( opts,'SCRF_HIG' ) || isempty( opts.SCRF_HIG ), opts.SCRF_HIG = 0.90;                    end
  % approximation options
  if ~isfield( opts,'DEG_MAX'  ) || isempty( opts.DEG_MAX  ), opts.DEG_MAX  = 128;                     end
  if ~isfield( opts,'m_r_arn'  ) || isempty( opts.m_r_arn  ), opts.m_r_arn  = min( 64, opts.DEG_MAX ); end
  if ~isfield( opts,'r_arn'    ) || isempty( opts.r_arn    ), opts.r_arn    = min( opts.m_r_arn, 64 ); end
  if ~isfield( opts,'r_lan'    ) || isempty( opts.r_lan    ), opts.r_lan    = 64;                      end
  % matrix options
  if ~isfield( opts,'skew'     ) || isempty( opts.skew     ), opts.skew     = 0; info.A.skew = 0;      else info.A.skew = opts.skew; end
  if ~isfield( opts,'nstrips'  ) || isempty( opts.nstrips  ), opts.nstrips  = 2;                       end
  if ~isfield( opts,'fov_tol'  ) || isempty( opts.fov_tol  ), opts.fov_tol  = 1e-3;                    end

  %% 03. PROBLEM SETUP AND MATRIX ANALYSIS: FIND FOV AND RITZ VALUES
  % augment matrix A
  f  = u( :,1 );  u( :,1 ) = [];
  f_ = flip( eye( p, 1 ) );
  idu = find( flip( any( u ) ) );
  u = u( :, flip( find( any( u ) ) ) );
  function [ x, x_ ] = Atilde( x, x_ )
    x = A( x ) + u * x_( idu );
    for i = 1 : p - 1;
      x_( i ) = x_( i + 1 );
    end
    x_( p ) = 0;
  end
  % (iom)-Arnoldi / Lanczos calls
  make_points  = ~isfield( info.A,'pts'     ) || isempty( info.A.pts     );
  make_polygon = ~isfield( info.A,'polygon' ) || isempty( info.A.polygon );
  info = bamphi_fov( n, A, At, @Atilde, f, p, make_points, make_polygon, info, opts );

  %% 04. EXECUTION

  %% Initialize history
  if ~isfield( info.A, 'his' ) || isempty( info.A.his )
    info.A.his.tau = []; info.A.his.p = []; info.A.his.tol = []; % initialize history
  end

  % Timesteps analysis
  T = [ 0; t(:) ];
  if ( length( T ) == 2 )
    predec = [ 0 1 ];
  else
    [ th,rh ] = cart2pol( real( t ), imag( t ) );
    if all( abs( mod( th, pi ) - mod( th( 1 ), pi ) ) < 1e-3 )
      % they're aligned in the complex plane
      idx_1 = find( abs( th - th( 1 ) ) <  1e-3 ); [ ~, id_1 ] = sort( rh( idx_1 ) );
      idx_2 = find( abs( th - th( 1 ) ) >= 1e-3 ); [ ~, id_2 ] = sort( rh( idx_2 ) );
      predec( idx_1( id_1 ) + 1 ) = [ 1, idx_1( id_1( 1:end-1 ) ) + 1 ];
      predec( idx_2( id_2 ) + 1 ) = [ 1, idx_2( id_2( 1:end-1 ) ) + 1 ];
    else
      % they're scattered in the complex plane
      [ ~, predec ] = minspantree( graph( abs( repmat( T,1,length( T ) ) - repmat( T,1,length( T ) ).' ) ) );
    end
  end
  done = [ true, false( 1, length( predec ) - 1 ) ];

  % Execution
  for j = 2 : length( T )
    [ info, f, f_, done ] = timesteps_manager( j, T, predec, @Atilde, A, u, idu, p, opts, info, f, f_, done );
  end

  %% 05. ADJUST AND RETURN
  f = n_u * f( :, 2 : length( T ) );

end

%
function [ info, f, f_, done ] = timesteps_manager( j, T, predec, Atilde, A, u, idu, p, opts, info, f, f_, done )

  if done( j )
    return
  end

  % Recursive bit
  % (I know you don't like it but it NEVER occurs if you don't ask for timesteps
  % scattered in the complex plane)
  if ~done( predec( j ) )
    [ info, f, f_, done ] = timesteps_manager( predec( j ), T, predec, Atilde, A, u, idu, p, opts, info, f, f_, done );
    done( predec( j ) ) = true;
  end

  % Perform Backward Error Analysis
  tau = T( j ) - T( predec( j ) );
  id = find( ( info.A.his.tau == tau ) .* ( info.A.his.p == p ) .* ( info.A.his.tol == opts.tol ), 1, 'first' );
  if isempty( id )

    id = length( info.A.his.tau ) + 1;
    info.A.his.tau( id ) = tau;
    info.A.his.p  ( id ) = p;
    info.A.his.tol( id ) = opts.tol;
    info.A.his.ell( id ) = max( opts.DEG_MAX - info.A.his.p( id ) - info.A.pts.r, 1 );
    if p
      info.A.his.x{ id } = [ zeros( 1, p ), info.A.pts.rho, info.A.pts.mu * ones( 1,info.A.his.ell( id ) ) ];
    else
      info.A.his.x{ id } = [ info.A.pts.mu, info.A.pts.rho, info.A.pts.mu * ones( 1,info.A.his.ell( id ) - 1 ) ]; % gotta have a starter for talezer
    end

    info = bamphi_bea( id, tau, info, opts );

    info.A.his.map{ id } = 2 * ( imag( info.A.his.x{ id }( : ) ) ~= 0 ) - 1;
    info.A.his.map{ id } = [ info.A.his.map{ id }( 1 ); info.A.his.map{ id }; -info.A.his.map{ id }( end ) ];

    % Scaling Refinement Params
    info.A.his.its_s{ id } = zeros( info.A.his.s( id ), 1 );
    info.A.his.reject( id ) = 0;

  end

  % PERFORM KRYLOV/NEWTON-HERMITE APPROXIMATION
  if not( isfield( info,'scal_ref' ) )
    id_scal_ref = 1;
    info.scal_ref.ts  ( id_scal_ref ) = info.A.his.tau( id ) / info.A.his.s( id );
    info.scal_ref.tol ( id_scal_ref ) = info.A.his.tol( id );
    info.scal_ref.jump( id_scal_ref ) = 1;
  else
    id_scal_ref = find( ( info.scal_ref.tol == info.A.his.tol( id ) ), 1, 'first' );
    if isempty( id_scal_ref )
      id_scal_ref = length( info.scal_ref.ts ) + 1;
      info.scal_ref.tol ( id_scal_ref ) = info.A.his.tol( id ); % for this tolerance...
      info.scal_ref.ts  ( id_scal_ref ) = info.A.his.tau( id ) / info.A.his.s( id ); % for this measure unit...
      info.scal_ref.jump( id_scal_ref ) = 1; % ...we jump this way
    end
  end
  jump = max( 1, round( info.scal_ref.jump( id_scal_ref ) * info.scal_ref.ts( id_scal_ref ) * info.A.his.s( id ) / info.A.his.tau( id ) ) );
  info.scal_ref.ts  ( id_scal_ref ) = info.A.his.tau( id ) / info.A.his.s( id ); % for this measure unit...
  info.scal_ref.jump( id_scal_ref ) = jump; % ...we jump this way
  step = 1;
  % Krylov step
  jump_lim = inf;
  SCRF_LOW = opts.SCRF_LOW;
  if not( isfield( info.A, 'arn' ) && isfield( info.A.arn, 'done' ) && not( info.A.arn.done ) )
     w =  f( :, predec( j ) );
    w_ = f_( :, predec( j ) );
  else
    while isfield( info.A, 'arn' ) && isfield( info.A.arn, 'done' ) && not( info.A.arn.done ) %&& false

      info.scal_ref.jump( id_scal_ref ) = jump;
      if info.A.his.s( id ) - step + 1 < jump
        jump = info.A.his.s( id ) - step + 1;
      end
      if not( opts.scal_ref )
        jump = 1;
      end
      info.A.his.ts( id ) = jump * info.A.his.tau( id ) / info.A.his.s( id );
      l = length( info.A.his.d{ id } );
      while l < jump
        info.A.his.d{ id }{ l+1 } =                          info.A.his.d{ id }{ l   }   .* [ 1; cumprod(   l     * ones( length( info.A.his.d{ id }{ l   } ) - 1 ,1 ) ) ];
        info.A.his.d{ id }{ l+1 } = ( info.A.his.F{ id }.' * info.A.his.d{ id }{ l+1 } ) ./ [ 1; cumprod( ( l+1 ) * ones( length( info.A.his.d{ id }{ l+1 } ) - 1 ,1 ) ) ];
        l = l + 1;
        if any( isinf( info.A.his.d { id }{ l } ) ) || any( isnan( info.A.his.d { id }{ l } ) )
          jump     = l - 1;
          jump_lim = l - 1;
          info.A.his.d { id }( l ) = [];
          break;
        end
      end

      rly = ( info.A.pts.rly && isreal( f ) && isreal( info.A.his.ts( id ) ) );
      [ w, w_, its, e_t, ret ] = bamphi_krylovich( 'assemble', jump, [], A, [], p, id, rly, opts, info ); info.A.arn.done = true;
      info.A.arn.e_t = e_t;
      info.A.his.its_s{ id }( step ) = its + 1;

      step = step + jump;

      if opts.scal_ref && ( info.A.his.s( id ) > 1 )
        if not( ret )
          if ( e_t < SCRF_LOW * opts.DEG_MAX )
            jump = jump + 1;
          end
          if ( e_t > opts.SCRF_HIG * opts.DEG_MAX ) && ( jump > 1 )
            jump = jump - 1;
          end
          SCRF_LOW = min( opts.SCRF_LOW, 1.05 * SCRF_LOW );
        end
        if ret && ( jump > 1 )
          step = step - jump;
          jump = jump - 1;
          info.A.his.reject( id ) = info.A.his.reject( id ) + 1;
          SCRF_LOW = 0.5 * SCRF_LOW;
          info.A.arn.done = false;
        end
      end
      jump = min( jump, jump_lim );
    end% endwhile
  end % endif

  % Hermite-Newton step(s)
  while step < info.A.his.s( id ) + 1 % true execution

    if opts.scal_ref
      v  = w;
      v_ = w_;
    end

    % begin only important if scal_ref is on
    info.scal_ref.jump( id_scal_ref ) = jump;
    if info.A.his.s( id ) - step + 1 < jump
      jump = info.A.his.s( id ) - step + 1;
    end
    if not( opts.scal_ref )
      jump = 1;
    end
    info.A.his.ts( id ) = jump * info.A.his.tau( id ) / info.A.his.s( id );
    l = length( info.A.his.d { id } );
    while l < jump
      info.A.his.d { id }{ l+1 } =                          info.A.his.d{ id }{ l   }   .* [ 1; cumprod(   l     * ones( length( info.A.his.d{ id }{ l   } ) - 1 ,1 ) ) ];
      info.A.his.d { id }{ l+1 } = ( info.A.his.F{ id }.' * info.A.his.d{ id }{ l+1 } ) ./ [ 1; cumprod( ( l+1 ) * ones( length( info.A.his.d{ id }{ l+1 } ) - 1 ,1 ) ) ];
      l = l + 1;
      if any( isinf( info.A.his.d { id }{ l } ) ) || any( isnan( info.A.his.d { id }{ l } ) )
        jump     = l - 1;
        jump_lim = l - 1;
        info.A.his.d { id }( l ) = [];
        break;
      end
    end
    % end only important if scal_ref is on

    rly = ( info.A.pts.rly && isreal( f ) && isreal( info.A.his.ts( id ) ) );
    [ w, w_, its, ret ] = bamphi_talezervich( id, jump, opts, info, Atilde, A, w, w_, rly, 0 ); % 2 / 3 * SCRF_LOW * opts.DEG_MAX
    info.A.his.its_s{ id }( step ) = its + 1;
    info.A.his.e_t = its + 1;

    step = step + jump;

    if opts.scal_ref && ( info.A.his.s( id ) > 1 )
      if not( ret )
        if ( its < SCRF_LOW * opts.DEG_MAX )
          jump = jump + 1;
        end
        if ( its > opts.SCRF_HIG * opts.DEG_MAX ) && ( jump > 1 )
          jump = jump - 1;
        end
        SCRF_LOW = min( opts.SCRF_LOW, 1.05 * SCRF_LOW );
      end
      if ret && ( jump > 1 )
        step = step - jump;
        jump = jump - 1;
        info.A.his.reject( id ) = info.A.his.reject( id ) + 1;
        SCRF_LOW = 0.5 * SCRF_LOW;
         w = v;
        w_ = v_;
      end
    end
    jump = min( jump, jump_lim );

  end %endwhile

  info.scal_ref.SCRF_LOW = SCRF_LOW;

   f( :, j ) = w;
  f_( :, j ) = w_;
  done( j )  = true;

end

function out = my_norm( in, opts, flag )

  if opts.earl_ter && flag
    out = norm( in, opts.norm );
  else
    out = [];
  end

end

%
function [ p_A, p_A_, its, ret ] = bamphi_talezervich( id, jump, opts, info, Atilde, A, w1, w1_, rly, e_t )
%
%
  ret = 1;
  % Set initial parameters
  O_or_1i = 1i * ( 1 - rly );
  m = length( info.A.his.x{ id } ) - 1;
  cof_w1 = zeros( m + 1, 1 );
  cof_w2 = zeros( m + 1, 1 );
  c1 = 0;
  c2 = 0;
  tol = info.A.his.tol( id );
  ts = info.A.his.ts( id );

  % Separate real and imaginary part to work with reals only if possible
  z_r = real( info.A.his.x{ id } );         z_i = imag( info.A.his.x{ id } );
  d_r = real( info.A.his.d  { id }{ jump } ); d_i = imag( info.A.his.d  { id }{ jump } );
  cof_w1 = d_r + O_or_1i * d_i;

  % First p + 1 iterations
  k = 1;
  w1  = cof_w1( k ) * w1;
  w1_ = cof_w1( k ) * w1_;
  p_A  = w1 ;
  p_A_ = w1_;
  while sqrt( w1_' * w1_ )
    % Coefficient and Product
    [ w1, w1_ ] = Atilde( w1, w1_ );
    w1  = ( ts * cof_w1( k + 1 ) / cof_w1( k ) ) * w1;
    w1_ = ( ts * cof_w1( k + 1 ) / cof_w1( k ) ) * w1_;
    p_A  = p_A  + w1;
    p_A_ = p_A_ + w1_;
    k = k + 1;
  end
  c2 = my_norm( w1, opts, true );

  if not( rly ) || all( info.A.his.map{ id }(1:end-1) == -1 )
    % Successive Iterations
    % everything's real or complex and unpaired: go with Newton
    while ( k < m + 1 )
      % Coefficient and Product
      w1 = ( ts * cof_w1( k + 1 ) / cof_w1( k ) ) * ( A( w1 ) - info.A.his.x{ id }( k ) * w1 );
      c1 = c2;
      c2 = my_norm( w1, opts, k > e_t );
      if not( isempty( c2 ) ) && ( isinf( c2 ) || isnan( c2 ) )
        break, % emergency break
      end
      p_A = p_A + w1;
      k = k + 1;
      % Early termination
      if ( opts.error )
        % relative error
        if ( c1 + c2 <= tol * my_norm( p_A, opts, k > e_t ) )
          ret = 0;
          break
        end
      else
        % absolute error
        if ( c1 + c2 <= tol )
          ret = 0;
          break
        end
      end
    end%while

  else
    cof_w2 = d_r;
    % Successive Iterations
    % we got complex data points and they're in conjugated pairs: go with Tal-Ezer
    while ( k < m + 1 )
      if ( info.A.his.map{ id }( k + 1 ) == -1 ) || not( rly )
        % Coefficient and Product
        w1 = ( ts * cof_w1( k + 1 ) / cof_w1( k ) ) * ( A( w1 ) - info.A.his.x{ id }( k ) * w1 );
        c1 = c2;
        c2 = my_norm( w1, opts, k > e_t );
        if not( isempty( c2 ) ) && ( isinf( c2 ) || isnan( c2 ) )
          break, % emergency break
        end
        p_A = p_A + w1;
        k = k + 1;
      else
        % Coefficient and Product
        w2 = ( ts * cof_w2( k + 1 ) / cof_w1( k ) ) * ( A( w1 ) - z_r( k ) * w1 );
        if ( ( info.A.his.map{ id }( k + 1 ) * info.A.his.map{ id }( k ) ) == -1 )
          p_A = p_A + w2;
          c1 = c2;
        else
          p_A = p_A + w2 + w1;
          c1 = my_norm( w1, opts, k > e_t );
        end
        c2 = my_norm( w2, opts, k > e_t );
        % Coefficient and Product
        w1 = ( ts * cof_w1( k + 2 ) / cof_w2( k + 1 ) ) * ...
             ( A( w2 ) - z_r( k + 1 ) * w2 + ( ts * z_i( k + 1 )^2 * cof_w2( k + 1 ) / cof_w1( k ) ) * w1 );
        if ( k < m ) && ( ( info.A.his.map{ id }( k + 3 ) * info.A.his.map{ id }( k + 2 ) ) == -1 )
          c1 = c2;
          c2 = my_norm( w1, opts, k > e_t );
          if not( isempty( c2 ) ) && ( isinf( c2 ) || isnan( c2 ) )
            break, % emergency break
          end
          p_A = p_A + w1;
        end
        k = k + 2;
      end %endif
      if ( opts.error )
        % relative error
        if ( c1 + c2 <= tol * my_norm( p_A, opts, k > e_t ) )
          ret = 0;
          break
        end
      else
        % absolute error
        if ( c1 + c2 <= tol )
          ret = 0;
          break
        end
      end

    end % endwhile
  end
  its = k - 1;

end
%endfunction
