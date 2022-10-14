function info = bamphi_bea( id, tau, info, opts )
% BEA_CC_POLYCHAIN given the vertices of a convex compact set K s.t. 0 \in K and
% K contains W( A ), return the scaling parameter s and the divided differences
% d_s at x * ( t / s ).
%
% Author: Franco Zivcovich
% Last Update: January the 6th, 2022.
%

  n_sample = 64;

  x   = info.A.his.x{ id }(:) - info.A.pts.mu; % just for computing divided differences quickly later on
  p   = info.A.his.p( id );
  r   = info.A.pts.r;
  ell = info.A.his.ell( id );
  m   = r + p + ell;
  L   = 25;

  thresh = opts.tol / ( 1 + sqrt( 2 ) );

  lg2exp = log2( exp( 1 ) );
  lg2thr = log2( thresh );
  lg2fac = gammaln( m + 1 ) * lg2exp;

  % Set iteration
  s = 1; it = 0; max_it = 1e6;

  % Quick iteration
  s = s - 1; err = inf;

  z = cohupts( info.A.polygon, n_sample );

  lg2z = sum( log2( z * ones( 1,r ) - ones( n_sample,1 ) * info.A.pts.rho ), 2 ) + ( ell - 1 ) * log2( z - info.A.pts.mu ) + p * log2( z ); % this is log2( prod( x - x_i ) ) done cleverly
  nu = ( z * r - sum( info.A.pts.rho ) ) / ( m + 1 );
  while true

    while ( any( err > lg2thr ) || any( isnan( err ) ) || any( isinf( err ) ) ) && ( it < max_it )
      it = it + 1; s = s + 1; ts = info.A.his.tau( id ) / s;
      err = 2 * real( lg2exp * ts * nu( 1 ) - lg2fac + ( m - 1 ) * log2( ts ) + lg2z( 1 ) );
    end

    err = 2 * real( lg2exp * ts * nu - lg2fac + ( m - 1 ) * log2( ts ) + lg2z );

    if ( any( err > lg2thr ) || any( isnan( err ) ) || any( isinf( err ) ) ) && ( it < max_it )
      continue
    end
    break

  end

  % Thorough iteration
  s = s - 1; err = inf;
  lg2z = lg2z + log2( z - info.A.pts.mu );
  while true

    it = it + 1; s = s + 1; ts = info.A.his.tau( id ) / s;

    if opts.scal_ref
      [ d_s, F_s, s_d ] = bamphi_dd_phi( [ ts * x; zeros( L, 1 ) ], 0, L + ell - 1 );
    else
        d_s             = bamphi_dd_phi( [ ts * x; zeros( L, 1 ) ], 0, L + ell - 1 );
    end
    d_s = exp( ts * info.A.pts.mu ) * d_s;
    err = log2( exp( - ts * z ) .* ( bsxfun( @power, ts * z, 0 : L ) * d_s( m - 1 : m + L - 1 ) ) ) + m * log2( ts ) + lg2z;
    err = abs( log1p( - 2.^err ) ./ ( ts * ( z - info.A.pts.mu ) ) );

    if ( any( err > thresh ) || any( isnan( err ) ) || any( isinf( err ) ) ) && ( it < max_it )
      continue
    end
    break

  end
  if not( it < max_it ),
    error('bamphi_bea: unable to determine a scaling strategy');
  end

  info.A.his.s( id ) = s;
  info.A.his.d{ id }{ 1 } = d_s( 1 : m );
  if opts.scal_ref
    info.A.his.F{ id } = exp( ts * info.A.pts.mu ) * ( F_s( 1 : m, 1 : m )^s_d );
  end

end

%
function p = cohupts( k, N )
% Return N points from the border of the convex polygon t * polygon / s.
  R = randi( length( k ), [ N,1 ] ); r = rand( N,1 );
  p = r .* k( R ) + ( 1 - r ) .* k( mod( R,length( k ) ) + 1 );
end
