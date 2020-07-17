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
    
function [xs,ws] = gauss_legendre_rule(N,a,b)

mean = (a+b)/2; scale = (b-a)/2;

% find Gauss-Legendre rule for int(scale*f(scale*y+mean),y,-1,1), and set
% wi' = wi*scale, xi' = scale*xi+mean
alphas = zeros(1,N);
betas = (1:N)./sqrt((2.*(1:N)-1).*(2.*(1:N)+1));
J = sparse([(1:N),(1:N-1),(2:N)], [(1:N),(2:N),(1:N-1)], [alphas,betas(1:N-1),betas(1:N-1)], N, N);

xs = eig(J); % get eigenvalues of sparse matrix (fast)
betas = [1,betas];

% get eigenvectors (unnormalized)
bqs = zeros(length(xs),N);
bqs(:,1) = 1;
qs = zeros(length(xs),N);
qs(:,1) = 1;
if N > 1
    bqs(:,2) = (xs-alphas(1));
    qs(:,2) = bqs(:,2)./betas(2);
    for m = 3:N
        bqs(:,m) = (xs-alphas(m-1)).*bqs(:,m-1)./betas(m-1)-betas(m-1)*bqs(:,m-2)/betas(m-2);
        qs(:,m) = bqs(:,m)./betas(m);
    end
end
qnorm2s = sum(qs.^2,2);
xs = xs*scale+mean;
ws = (2*scale)./qnorm2s;
ws(isnan(ws)) = 0;

end
