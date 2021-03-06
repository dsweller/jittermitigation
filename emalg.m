% Copyright (c) 2011, Daniel Weller
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

% This code is based on the EM algorithm code from:
% 
% Weller, Daniel S. and Goyal, Vivek K. "On the estimation of nonrandom
%      signal coefficients from jittered samples." IEEE Trans. Signal
%      Process., vol. 59, no. 2, pp. 587-597, Feb. 2011. DOI:
%      10.1109/TSP.2010.2090347

function [x_hat, loglikelihood, diff, ldiff, iter] = emalg(y, x_hat_0, sigma_ws, weights_sigma_w, weights_sigma_z, weights_z, Hzs)

K = length(x_hat_0);
N = length(y);
J1 = length(weights_sigma_w);
J2 = length(weights_sigma_z);
J3 = size(weights_z,2);

% iterate using EM algorithm
maxiters = 100;
delta = 1e-8;
epsilon = 1e-8;

% when computing flops, ignore rescaling used to mitigate overflow/underflow

y_minus_Hzxs = repmat(y,[1,J2,J3]);
for k = 1:K
    y_minus_Hzxs = y_minus_Hzxs - Hzs(:,:,:,k).*x_hat_0(k);
end

logs = y_minus_Hzxs.^2;

logs_scalefactor = min(min(logs./max(2.*sigma_ws.^2) - repmat(shiftdim(log(weights_z),-1),[N,1,1]),[],3),[],2);
logs_scalefactor_rep = repmat(logs_scalefactor,[1,J2,J3]);
exps = 0;
for j1 = 1:J1
    exps = exps + (weights_sigma_w(j1)/sqrt(2*pi*sigma_ws(j1)^2)).*exp(logs_scalefactor_rep-logs./(2*sigma_ws(j1)^2));
end

likelihoods = 0;
for j2 = 1:J2
    likelihoods = likelihoods + weights_sigma_z(j2).*sum(exps(:,j2,:).*repmat(shiftdim(weights_z(j2,:),-1),[N,1,1]),3);
end

loglikelihood = sum(log(likelihoods)-logs_scalefactor);

x_hat = x_hat_0;

diff = Inf;
ldiff = Inf;

Aadd = zeros([N,K,K]);
for iter = 1:maxiters
    % compute E[H(z)|Y;x] and E[H(z)'H(z)|Y;x]
    exps = 0;
    for j1 = 1:J1
        exps = exps + (weights_sigma_w(j1)/sqrt(2*pi*sigma_ws(j1)^2)/(sigma_ws(j1)^2)).*exp(logs_scalefactor_rep-logs./(2*sigma_ws(j1)^2));
    end
    
    F = 0;
    A = 0;
    for j2 = 1:J2
        exps_weighted = reshape(exps(:,j2,:),[N,J3]).*weights_z(j2+zeros(N,1),:);
        Hzs_j2 = reshape(Hzs(:,j2,:,:),[N,J3,K]);
        exps_weighted_Hs = exps_weighted(:,:,ones(K,1)).*Hzs_j2;
        F = F + weights_sigma_z(j2).*sum(exps_weighted_Hs,2);
        for k2=1:K
            Hzs_j2_k = Hzs_j2(:,:,k2);
            for k1=1:K
                Aadd(:,k1,k2) = sum(exps_weighted_Hs(:,:,k1).*Hzs_j2_k,2);
            end
        end
        A = A + weights_sigma_z(j2).*Aadd;
    end
    
    F = reshape(F,[N,K])./repmat(likelihoods,[1,K]);
    A = reshape(A,[N,K,K])./repmat(likelihoods,[1,K,K]);
    A = reshape(sum(A,1),K,K);
    
    B = (F'*y);

    % compute new x
    old_x_hat = x_hat;
    x_hat = (A+1e-9.*eye(K))\B;
    
    % compute new log-likelihood
    y_minus_Hzxs = repmat(y,[1,J2,J3]);
    for k = 1:K
        y_minus_Hzxs = y_minus_Hzxs - Hzs(:,:,:,k).*x_hat(k);
    end
    
    logs = y_minus_Hzxs.^2;

    % again, ignore rescaling when computing op count
    logs_scalefactor = min(min(logs./max(2.*sigma_ws.^2) - repmat(shiftdim(log(weights_z),-1),[N,1,1]),[],3),[],2);
    logs_scalefactor_rep = repmat(logs_scalefactor,[1,J2,J3]);
    exps = 0;
    for j1 = 1:J1
        exps = exps + (weights_sigma_w(j1)/sqrt(2*pi*sigma_ws(j1)^2)).*exp(logs_scalefactor_rep-logs./(2*sigma_ws(j1)^2));
    end
    
    likelihoods = 0;
    for j2 = 1:J2
        likelihoods = likelihoods + weights_sigma_z(j2).*sum(exps(:,j2,:).*repmat(shiftdim(weights_z(j2,:),-1),[N,1,1]),3);
    end
    
    old_loglikelihood = loglikelihood;
    loglikelihood = sum(log(likelihoods)-logs_scalefactor);

    if loglikelihood < old_loglikelihood % maximize-step failed; exit
        loglikelihood = old_loglikelihood;
        x_hat = old_x_hat;
        break;
    end
    
    diff = sum((x_hat - old_x_hat).^2);
    ldiff = loglikelihood - old_loglikelihood;
    
    if diff < delta || ldiff < epsilon % convergence; exit
        break;
    end
end
% disp(sprintf('Difference = %g (%d iterations)', diff, iter));

end
