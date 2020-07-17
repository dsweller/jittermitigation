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

% weighting function:
% w(x) = (1/sqrt(2*pi*sigma^2))exp(-(x-mu)^2/(2*sigma^2))
%
% int(f(x)w(x),x,-Inf,Inf) =~ sum(ws(i)*f(xs(i)))
%
function [xs, ws] = gauss_hermite_rule(N, mu, sigma)

betas = sqrt((1:N-1));
J = sparse([(1:N-1),(2:N)], [(2:N),(1:N-1)], [betas(1:N-1),betas(1:N-1)], N, N);

if (N > 500)
    [V,D] = eig(full(J));
    xs = diag(D)*sigma+mu;
    ws = V(1,:)'.^2;
else
    xs = eig(J); % get eigenvalues of sparse matrix (fast)

    % get eigenvectors (unnormalized)
    qs = zeros(length(xs),N);
    qs(:,1) = 1;
    if N > 1
        qs(:,2) = xs;
        for m = 3:N
            qs(:,m) = (xs.*qs(:,m-1)-betas(m-2)*qs(:,m-2))/betas(m-1);
        end
    end
    qnorm2s = sum(qs.^2,2);
    xs = xs*sigma+mu;
    ws = 1./qnorm2s;
end

end
