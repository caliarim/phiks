function varargout = wantcache(varargin)
% WANTCACHE -- Return or update Expint cache desirabilty status
%
% SYNOPSIS:
%       wantcache(state)
%   b = wantcache
%
% PARAMETERS:
%   state - New cache desirability state.  OPTIONAL.
%           Used as a boolean flag stating whether or not caching is
%           wanted in PHIPADE and ONEEZ2.
%
%           Possible values are
%             - FALSE, 0, 'false', 'off', 'no', 'n' -- caching unwanted
%             - anything else                       -- caching wanted
%
%           Case is insignificant in string values.
%
% RETURNS:
%   b   - Boolean value indicating whether or not caching is wanted.
%         OPTIONAL.
%
% SEE ALSO:
%   TRUE, FALSE, PHIPADE, ONEEZ2.

% This file is part of the 'Expint'-package,
% see http://www.math.ntnu.no/num/expint/
%
% $Revision: 1.5 $  $Date: 2005/11/09 16:37:20 $

error(nargchk(0, 1, nargin));
error(nargoutchk(0, 1, nargout));

if all([nargin, nargout] == 0),
   error('wantcache:einval', 'Incorrect usage, see documentation');
end

persistent state;

% Default state is to turn caching on:
if isempty(state),
   state = true; 
end

% Set new state:
if nargin > 0,
   if islogical(varargin{1}),
      state = varargin{1};
   else
      switch lower(varargin{1}),
      case { 0, 'false', 'off', 'no', 'n' },
         state = false;
      otherwise
         state = true;
      end
   end
end

% Possibly return current state:
if nargout > 0,
   varargout{1} = state;
end
