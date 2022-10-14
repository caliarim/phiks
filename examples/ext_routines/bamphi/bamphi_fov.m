function info = bamphi_fov( n, A, At, Atilde, f, p, make_points, make_polygon, info, opts )
%
%
%
%
  if ( nargin < 10 ),                            opts = [];                  end
  if ( nargin <  9 ),                            info = [];                  end
  if ( nargin <  8 || isempty( make_polygon ) ), make_polygon = true;        end
  if ( nargin <  7 || isempty( make_points  ) ), make_points  = true;        end
  if ( nargin <  6 || isempty( p ) ),            p = 0;                      end
  if ( nargin <  5 || isempty( f ) ),            f = randn( n, 1 );          end
  if ( nargin <  4 || isempty( Atilde ) ),       Atilde = @( x, x_ ) A( x ); end
  if ( nargin <  3 ),                            At = [];                    end
  if ( nargin <  2 ), error('bamphi_fov: at least two inputs are required'); end

  % % technical options
  if ~isfield( opts,'low_stor' ) || isempty( opts.low_stor ), opts.low_stor = true;                    end
  if ~isfield( opts,'iom'      ) || isempty( opts.iom      ), opts.iom      = true;                    end
  % % approximation options
  if ~isfield( opts,'DEG_MAX'  ) || isempty( opts.DEG_MAX  ), opts.DEG_MAX  = 128;                     end
  if ~isfield( opts,'m_r_arn'  ) || isempty( opts.m_r_arn  ), opts.m_r_arn  = min( 64, opts.DEG_MAX ); end
  if ~isfield( opts,'r_arn'    ) || isempty( opts.r_arn    ), opts.r_arn    = min( opts.m_r_arn, 64 ); end
  if ~isfield( opts,'k_lan'    ) || isempty( opts.k_lan    ), opts.k_lan    = 64;                      end
  % % matrix options
  if ~isfield( opts,'skew'     ) || isempty( opts.skew     ), opts.skew     = 0; info.A.skew = 0;      else info.A.skew = opts.skew; end
  if ~isfield( opts,'nstrips'  ) || isempty( opts.nstrips  ), opts.nstrips  = 2;                       end
  if ~isfield( opts,'fov_tol'  ) || isempty( opts.fov_tol  ), opts.fov_tol  = 1e-3;                    end

  T_r = []; T_i = [];
  if ( make_points || make_polygon )
    fov = zeros( 2 );
    if ( make_polygon && not( isempty( At ) ) )
      v = randn( n, 1 );
      if ( info.A.skew ~= -1 ), [ T_r, ~, fov( 1,: ), it, v ] = lanczos_fov( 1/2,  1/2, A, At, v, opts.k_lan,  1, opts.fov_tol, opts.low_stor ); end % horizontal strip
      if ( info.A.skew ~=  1 ), [ T_i, ~, fov( 2,: ), it, v ] = lanczos_fov( 1/2, -1/2, A, At, v, opts.k_lan, -1, opts.fov_tol, opts.low_stor ); end % vertical   strip
      if not( info.A.skew ),
        info.A.skew = ( fov( 2,1 ) == fov( 2,2 ) ) - ( fov( 1,1 ) == fov( 1,2 ) );
      end
      clear v;
    end

    if ( make_points || isempty( At ) )

      if ( info.A.skew == 0 || ( isempty( T_r ) && isempty( T_i ) ) )
        H = bamphi_krylovich( 'build', [], Atilde, A, f, p, [], [], opts, info ); info.A.arn.done = opts.low_stor;
      elseif ( info.A.skew ==  1 )
        H = T_r; % A is hermitian
      elseif ( info.A.skew == -1 )
        H = T_i; % A is skew-hermitian
      end
      H = H( 1 : length( H ) - 1, 1 : length( H ) - 1 );

      % Build Hermite-Ritz interpolation sequence
      if make_points
        info.A.pts.rly = isreal( H );
        info.A.pts.rho = eig( H );
        if     ( norm( H - H',1 ) < max( abs( info.A.pts.rho ) ) * 1e-10 ) || ( info.A.skew ==  1 )
          info.A.pts.rho =      real( info.A.pts.rho );
        elseif ( norm( H + H',1 ) < max( abs( info.A.pts.rho ) ) * 1e-10 ) || ( info.A.skew == -1 )
          info.A.pts.rho = 1i * imag( info.A.pts.rho );
        end

        info.A.pts.mu  = mean( info.A.pts.rho );
        info.A.pts.rho = info.A.pts.rho - info.A.pts.mu;
        info.A.pts.rho = norm( info.A.pts.rho, inf ) * bamphi_hermite_ext( [], info.A.pts.rho( : ) / norm( info.A.pts.rho, inf ), info.A.pts.rly );
        info.A.pts.rho = info.A.pts.rho + info.A.pts.mu; % noshift
        info.A.pts.r   = length( info.A.pts.rho );
      end

    end

    if make_polygon && isempty( At )
      H_th = H + H'; rh = real( eig( H_th ) ); fov( 1,: ) = [ min( rh ), max( rh ) ] / 2; % horizontal strip
      H_th = H - H'; rh = imag( eig( H_th ) ); fov( 2,: ) = [ min( rh ), max( rh ) ] / 2; % vertical   strip
      if not( info.A.skew ),
        info.A.skew = ( fov( 2,1 ) == fov( 2,2 ) ) - ( fov( 1,1 ) == fov( 1,2 ) );
      end
    end

    % find rectangle's vertices
    info.A.polygon( 1,1 ) = fov( 1,1 ) + 1i * fov( 2,1 );
    info.A.polygon( 2,1 ) = fov( 1,2 ) + 1i * fov( 2,1 );
    info.A.polygon( 3,1 ) = fov( 1,2 ) + 1i * fov( 2,2 );
    info.A.polygon( 4,1 ) = fov( 1,1 ) + 1i * fov( 2,2 );

  end
end

%
function [ T, Q, mu, j, v ] = lanczos_fov( c, c_t, A, At, v, m, skew, tol, low_stor )
% Carries out m Complete Orthogonalization ARNOLDI iterations on the
% n - by - n matrix A and the vector q producing the Hessemberg matrix H.
%
% Author Franco Zivcovich 31 January 2021 franco.zivcovich@gmail.com
%

  j = 1;
  i_o = false;
  Q{ 1 } = v / sqrt( v' * v ); % v MUST be normalized outside
  while ( ( j < m + 1 ) && not( i_o ) )

    idj_m1 = ( mod( j - 2, 3 ) + 1 ) * low_stor + ( j - 1 ) * not( low_stor );
    idj    = ( mod( j - 1, 3 ) + 1 ) * low_stor +   j       * not( low_stor );
    idj_p1 = ( mod( j + 0, 3 ) + 1 ) * low_stor + ( j + 1 ) * not( low_stor );

    Q{ idj_p1 } = c * A( Q{ idj } ) + c_t * At( Q{ idj } );

    T( j,j ) = Q{ idj }' * Q{ idj_p1 };
    Q{ idj_p1 } = Q{ idj_p1 } - T( j,j ) * Q{ idj };
    if ( j > 1 )
      Q{ idj_p1 } = Q{ idj_p1 } - T( j-1,j ) * Q{ idj_m1 };
    end
    T( j+1,j ) = sqrt( Q{ idj_p1 }' * Q{ idj_p1 } );

    if ( j == 1 ) && ( T( j+1,j ) == 0 )
      mu = [ 0, 0 ];
      return
    end
    % check for convergence
    if ( j > 4 )
      [  S, mu ] = eig( T( 1:j,1:j ), 'vector' );
      [ mu, id ] = sort( ( skew == 1 ) * real( mu ) + ( skew == -1 ) * imag( mu ) );
      nm2_T = max( abs( mu ) );
      i_o = all( abs( S( j,[ id( 1 ), id( j ) ] ) ) * T( j+1,j ) / nm2_T < tol );
    end

    T( j,j+1 ) = skew * T( j+1,j );
    Q{ idj_p1 } = Q{ idj_p1 } / T( j+1,j );
    j = j + 1;
  end

  mu = [ mu( 1 ), mu( j-1 ) ];
  if not( low_stor )
    cf = sum( S( :,[ id( 1 ), id( j-1 ) ] ), 2 );
    v = cf( 1 ) * Q{ 1 };
    for i = 2 : j - 1
      v = v + cf( i ) * Q{ i };
    end
  end

end
%endfunction
