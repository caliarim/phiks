function suggestion = bamphi_suggest_arnoldi_size( opts, info )
% Given the info gathered, bamphi_arnoldi_size compute a possibly smart choice
% for next call's Krylov subspace size.

  if ~isfield( opts,'DEG_MAX'  ) || isempty( opts.DEG_MAX  ), opts.DEG_MAX  = 128;                     end
  if ~isfield( opts,'m_r_arn'  ) || isempty( opts.m_r_arn  ), opts.m_r_arn  = min( 64, opts.DEG_MAX ); end

  if isfield( info.A.his, 'e_t' )
    if isfield( opts,'r_arn' )
      aux = opts.r_arn;
    else
      aux = info.A.his.e_t;
    end
    suggestion = ceil( ( info.A.his.e_t + aux ) / 2 );
  else
    if isfield( opts,'r_arn' )
      aux = opts.r_arn;
    else
      aux = info.A.arn.e_t;
    end
    suggestion = ceil( ( info.A.arn.e_t + aux ) / 2 );
  end
  suggestion = min( [ opts.m_r_arn, suggestion ] );

end
