% Copyright (c) 2011, Daniel Weller (dweller@mit.edu)
% 
% If you use this code please cite:
% 
% Weller, Daniel S. and Goyal, Vivek K. "Bayesian post-processing methods
%      for jitter mitigation in sampling." IEEE Trans. Signal Process.,
%      vol. 59, no. 5, pp. 2112-2123, May 2011. DOI:
%      10.1109/TSP.2011.2108289
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%

% crossdiff.m -- helper function to compute differences between all x and y
% values; equivalent to repmat(x,1,length(y))-repmat(y,length(x),1) for
% vectors, but much faster
% 

function out = crossdiff(x,y)

out = x(:,ones(size(y))) - y(ones(size(x)),:);

end
