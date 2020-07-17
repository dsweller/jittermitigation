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
    
% This code performs performance comparisons between linear MMSE, ML, and
% Gibbs/slice sampler algorithms, such as those in Figure 10.

% This code uses the gammar() function from version 4.2 of the StatBox
% toolbox available for download at http://www.statsci.org/matlab/statbox.html
% or on the MATLAB Central Link Exchange. StatBox is copyright (c)
% 1996-2003, Gordon Smyth.
% 
% If the toolbox does not already reside on your MATLAB path, you need to
% adjust the below path to point to this toolbox.
addpath('../statbx42/');

% simulation parameters
numbootstrap = 500; % I_b for Algorithm 2 of paper
numiters = 500; % I for Algorithm 2 of paper
ntests = 1000; % # of trials for performance comparisons

% seed the randomizer (for both uniform and Gaussian generators)
randstate = 0;
randnstate = 0;
rand('twister',randstate); %#ok<RAND>
randn('state',randnstate); %#ok<RAND>

savepath = './performance.mat'; % path to final results

% signal parameters
K = 10; % # of parameters in signal model
h = @(t,K) sinc(t); % function handle for sinc basis evaluation
% mean_x_desired = 0;
sigma_x_desired = 1;
alpha_x = (K+3)/2; % see Eq. (10)
beta_x = (alpha_x-1).*sigma_x_desired.^2; % see Eq. (11)
shrinkagemethod = 'rej-midpt-th'; % hybrid rejection-midpoint-threshold method
initstrings = {'linear'}; % initialize with the no-jitter linear estimator

% observation model parameters
Ms = 4; % oversampling factor(s) for simulations
sigma_zs_desired = [0.5,0.25,0.2,0.15,0.1,0.05,0.025,0.02,0.015,0.01]; % desired (expected) jitter standard deviation
sigma_ws_desired = 0.05; % desired (expected) additive noise standard deviation

% other parameters
J1 = 9; J2 = 9; J3 = 129; % number of terms in quadratures (see Section II.A)

% temporary storage
sigma_xs_temp = zeros(ntests,1);
sigma_zs_temp = zeros(ntests,1);
sigma_ws_temp = zeros(ntests,1);
xs_temp = zeros(ntests,K);
x_ests_linunbiaseds_temp = zeros(ntests,K);
x_ests_linunbiased_knowns_temp = zeros(ntests,K);
x_ests_linnjs_temp = zeros(ntests,K);
x_ests_linnj_knowns_temp = zeros(ntests,K);
x_ests_ems_temp = zeros(ntests,K);
x_ests_em_knowns_temp = zeros(ntests,K);
x_ests_blsslices_temp = zeros(ntests,K);

% perform simulations
nsims = max([length(Ms),length(sigma_zs_desired),length(sigma_ws_desired)]);

% set parameters for all simulations
if length(Ms) < nsims, Ms = Ms(1).*ones(1,nsims); end
if length(sigma_zs_desired) < nsims
    sigma_zs_desired = sigma_zs_desired(1).*ones(1,nsims);
end
if length(sigma_ws_desired) < nsims
    sigma_ws_desired = sigma_ws_desired(1).*ones(1,nsims);
end
alpha_zs = (Ms.*K+3)./2;
beta_zs = (alpha_zs-1).*sigma_zs_desired.^2;
alpha_ws = (Ms.*K+3)./2;
beta_ws = (alpha_ws-1).*sigma_ws_desired.^2;

% permanent storage
sigma_xs = zeros(nsims,ntests,1);
sigma_zs = zeros(nsims,ntests,1);
sigma_ws = zeros(nsims,ntests,1);
xs = zeros(nsims,ntests,K);
x_ests_linunbiaseds = zeros(nsims,ntests,K);
x_ests_linunbiased_knowns = zeros(nsims,ntests,K);
x_ests_linnjs = zeros(nsims,ntests,K);
x_ests_linnj_knowns = zeros(nsims,ntests,K);
x_ests_ems = zeros(nsims,ntests,K);
x_ests_em_knowns = zeros(nsims,ntests,K);
x_ests_blsslices = zeros(nsims,ntests,K);

for is = 1:nsims
    M = Ms(is);
    alpha_z = alpha_zs(is);
    beta_z = beta_zs(is);
    sigma_z_desired = sigma_zs_desired(is);
    alpha_w = alpha_ws(is);
    beta_w = beta_ws(is);
    sigma_w_desired = sigma_ws_desired(is);

    N = M*K; % number of samples scales with oversampling factor and parameters

    % generate quadrature abscissas and weights (for use later)
    [abscissas_sigma_w,weights_sigma_w] = gen_laguerre_rule(J1,alpha_w-1,0,1);
    abscissas_sigma_w = abscissas_sigma_w(:).'./beta_w;
    weights_sigma_w = weights_sigma_w(:).'./sum(weights_sigma_w); % normalize
    [abscissas_sigma_z,weights_sigma_z] = gen_laguerre_rule(J2,alpha_z-1,0,1);
    abscissas_sigma_z = abscissas_sigma_z(:).'./beta_z;
    weights_sigma_z = weights_sigma_z(:).'./sum(weights_sigma_z); % normalize
    if sigma_z_desired > 0.1
        [abscissas_z,weights_z] = gauss_legendre_rule(J3,-pi/2,pi/2);
    else
        [abscissas_z,weights_z] = gauss_hermite_rule(J3,0,1);
    end

    % generate full set of Hzs (for use later)
    abscissas_zfull = zeros(J2,J3);
    weights_zfull = zeros(J2,J3);
    Hzs = zeros(N,J2,J3,K);
    for j2 = 1:J2
        sigma_z_j2 = sqrt(1./abscissas_sigma_z(j2));
        if sigma_z_desired > 0.1
            [abscissas_zfull(j2,:),weights_zfull(j2,:)] = tan_rule_transform(abscissas_z,weights_z,0,sigma_z_j2);
        else
            abscissas_zfull(j2,:) = abscissas_z.*sigma_z_j2;
            weights_zfull(j2,:) = weights_z;
        end
        Hzs(:,j2,:,:) = h(repmat(repmat((0:N-1).'./M,[1,1,J3])+repmat(shiftdim(abscissas_zfull(j2,:),-1),[N,1,1]),[1,1,1,K])-repmat(shiftdim(0:K-1,-2),[N,1,J3,1]),K);
    end

    % compute E[H(z)]
    EHz = 0;
    for j2 = 1:J2
        EHz = EHz + weights_sigma_z(j2).*sum(Hzs(:,j2,:,:).*repmat(shiftdim(weights_zfull(j2,:),-1),[N,1,1,K]),3);
    end
    EHz = reshape(EHz,[N,K]);

    % compute E[H(z)'H(z)]
    EHzTHz = 0;
    for n = 1:N
        Hzs_n = Hzs(n,:,:,:);
        Hzs_n_prod = permute(Hzs_n,[1,2,3,5,4]);
        Hzs_n_prod = Hzs_n(:,:,:,:,ones(K,1)).*Hzs_n_prod(:,:,:,ones(K,1),:);
        EHzTHz = EHzTHz + sum(repmat(weights_sigma_z,[1,1,1,K,K]).*sum(Hzs_n_prod.*repmat(shiftdim(weights_zfull,-1),[1,1,1,K,K]),3),2);
    end
    EHzTHz = reshape(EHzTHz,[K,K]);

    % compute E[H(z)H(z)']
    EHzHzT = EHz*EHz';
    diaginds = sub2ind([N,N],1:N,1:N);
    EHzHzT(diaginds) = sum(repmat(weights_sigma_z,[N,1]).*sum(sum(Hzs.^2,4).*repmat(shiftdim(weights_zfull,-1),[N,1,1]),3),2);

    % other initialization
    H0 = h(repmat((0:N-1).'./M,[1,K])-repmat(0:K-1,[N,1]),K);

    fprintf(1,'M = %d, desired sigma_z = %g, desired sigma_w = %g...', M, sigma_z_desired, sigma_w_desired); drawnow;

    % more temporary storage
    ys_temp = zeros(ntests,N);
    zs_temp = zeros(ntests,N);
    ws_temp = zeros(ntests,N);

    % generate LMMSE coefficients
    G_J = EHzHzT+(beta_w*(alpha_x-1)/beta_x/(alpha_w-1)).*eye(N);
    G_NJ = H0*H0'+(beta_w*(alpha_x-1)/beta_x/(alpha_w-1)).*eye(N);

    % do trials for this simulation
    for it = 1:ntests
        fprintf(1,'%4d/%4d',it,ntests); drawnow;

        % generate signal and samples
        sigma_x_temp = sqrt(beta_x./gammar(alpha_x));
        sigma_z_temp = sqrt(beta_z./gammar(alpha_z));
        sigma_w_temp = sqrt(beta_w./gammar(alpha_w));
        sigma_xs_temp(it) = sigma_x_temp;
        sigma_zs_temp(it) = sigma_z_temp;
        sigma_ws_temp(it) = sigma_w_temp;
        [x,y,z,w] = siggen(K,M,N,h,0,sigma_x_temp,sigma_z_temp,sigma_w_temp);
        xs_temp(it,:) = x;
        ys_temp(it,:) = y;
        zs_temp(it,:) = z;
        ws_temp(it,:) = w;    

        % now, take care of known Hzs
        if sigma_z_desired > 0.1 % use Legendre rule + transform
            [abscissas_zknown,weights_zknown] = gauss_legendre_rule(J3,-pi/2,pi/2);
            [abscissas_zknown,weights_zknown] = tan_rule_transform(abscissas_zknown,weights_zknown,0,sigma_z_temp);
            abscissas_zknown = abscissas_zknown(:).';
            weights_zknown = weights_zknown(:).';
        else % use Hermite rule
            [abscissas_zknown,weights_zknown] = gauss_hermite_rule(J3,0,sigma_z_temp);
            abscissas_zknown = abscissas_zknown(:).';
            weights_zknown = weights_zknown(:).';
        end
        Hzs_known = h(repmat(repmat((0:N-1).'./M,[1,1,J3])+repmat(shiftdim(abscissas_zknown,-1),[N,1,1]),[1,1,1,K])-repmat(shiftdim(0:K-1,-2),[N,1,J3,1]),K);

        % compute E[H(z)] for sigma_z known
        EHz_known = sum(Hzs_known.*repmat(shiftdim(weights_zknown,-1),[N,1,1,K]),3);
        EHz_known = reshape(EHz_known,[N,K]);

        % compute E[H(z)H(z)'] for sigma_z known
        EHzHzT_known = EHz_known*EHz_known';
        diaginds = sub2ind([N,N],1:N,1:N);
        EHzHzT_known(diaginds) = sum(sum(Hzs_known.^2,4).*repmat(shiftdim(weights_zknown,-1),[N,1,1]),3);

        % compute E[H(z)'H(z)] for sigma_z known
        EHzTHz_known = 0;
        for n = 1:N
            Hzs_n = Hzs_known(n,:,:,:);
            Hzs_n_prod = permute(Hzs_n,[1,2,3,5,4]);
            Hzs_n_prod = Hzs_n(:,:,:,:,ones(K,1)).*Hzs_n_prod(:,:,:,ones(K,1),:);
            EHzTHz_known = EHzTHz_known + sum(Hzs_n_prod.*repmat(shiftdim(weights_zknown,-1),[1,1,1,K,K]),3);
        end
        EHzTHz_known = reshape(EHzTHz_known,[K,K]);

        G_J_known = EHzHzT_known+(sigma_w_temp^2/sigma_x_temp^2).*eye(N);
        G_NJ_known = H0*H0'+(sigma_w_temp^2/sigma_x_temp^2).*eye(N);

        % generate linear estimators
        x_ests_linunbiaseds_temp(it,:) = EHz'*(G_J\y);
        x_ests_linunbiased_knowns_temp(it,:) = EHz_known'*(G_J_known\y);
        x_ests_linnjs_temp(it,:) = H0'*(G_NJ\y);
        x_ests_linnj_knowns_temp(it,:) = H0'*(G_NJ_known\y);

        % perform EM algorithm with unknown parameters
        x0 = H0\y; % use no-jitter initial value (classical linear est.)
        [x_hat,loglikelihood,diff,ldiff,iter] = emalg(y, x0(:), sqrt(1./abscissas_sigma_w), weights_sigma_w, weights_sigma_z, weights_zfull, Hzs);
        x_ests_ems_temp(it,:) = x_hat;

        % perform EM algorithm with known parameters
        x0 = H0\y; % use no-jitter initial value (classical linear est.)
        [x_hat,loglikelihood,diff,ldiff,iter] = emalg(y, x0(:), sigma_w_temp, 1, 1, weights_zknown, Hzs_known);
        x_ests_em_knowns_temp(it,:) = x_hat;

        % generate Gibbs-slice sampling Bayes MMSE estimate
        z0_temp = zeros(N,1);
        x0_temp = x_ests_linnjs_temp(it,:).'; % use no-jitter optimal linear estimate
        [z_est, x_est, sigma_x_est, sigma_w_est, sigma_z_est] = gibbsslicesampler(K, M, y, h, [], [], [], 'z0', z0_temp, 'x0', x0_temp, 'numiters', numiters, 'numbootstrap', numbootstrap, 'randomsigmax', {alpha_x,beta_x}, 'randomsigmaw', {alpha_w,beta_w}, 'randomsigmaz', {alpha_z,beta_z}, 'shrinkagemethod', shrinkagemethod);
        x_ests_blsslices_temp(it,:) = permute(x_est,[3,2,1]);
        fprintf(1,repmat('\b',1,9));
    end
    fprintf(1,'done\n'); drawnow;

    % store temporary results for permanent storage
    sigma_xs(is,:) = sigma_xs_temp;
    sigma_zs(is,:) = sigma_zs_temp;
    sigma_ws(is,:) = sigma_ws_temp;
    xs(is,:,:) = xs_temp;
    x_ests_linunbiaseds(is,:,:) = x_ests_linunbiaseds_temp;
    x_ests_linunbiased_knowns(is,:,:) = x_ests_linunbiased_knowns_temp;
    x_ests_linnjs(is,:,:) = x_ests_linnjs_temp;
    x_ests_linnj_knowns(is,:,:) = x_ests_linnj_knowns_temp;
    x_ests_ems(is,:,:) = x_ests_ems_temp;
    x_ests_em_knowns(is,:,:) = x_ests_em_knowns_temp;
    x_ests_blsslices(is,:,:,:) = x_ests_blsslices_temp;
end

% compute errors and approximate MSEs
e_linunbiaseds = x_ests_linunbiaseds - xs;
e_linunbiased_knowns = x_ests_linunbiased_knowns - xs;
e_linnjs = x_ests_linnjs - xs;
e_linnj_knowns = x_ests_linnj_knowns - xs;
e_ems = x_ests_ems - xs;
e_em_knowns = x_ests_em_knowns - xs;
e_blsslices = x_ests_blsslices - xs;
sqderr_linunbiaseds = sum(e_linunbiaseds.^2,3);
sqderr_linunbiased_knowns = sum(e_linunbiased_knowns.^2,3);
sqderr_linnjs = sum(e_linnjs.^2,3);
sqderr_linnj_knowns = sum(e_linnj_knowns.^2,3);
sqderr_ems = sum(e_ems.^2,3);
sqderr_em_knowns = sum(e_em_knowns.^2,3);
sqderr_blsslices = sum(e_blsslices.^2,3);
mse_linunbiaseds = mean(sqderr_linunbiaseds,2);
mse_linunbiased_knowns = mean(sqderr_linunbiased_knowns,2);
mse_linnjs = mean(sqderr_linnjs,2);
mse_linnj_knowns = mean(sqderr_linnj_knowns,2);
mse_ems = mean(sqderr_ems,2);
mse_em_knowns = mean(sqderr_em_knowns,2);
mse_blsslices = mean(sqderr_blsslices,2);

% compute confidence intervals
varse_linunbiaseds = var(sqderr_linunbiaseds,0,2);
varse_linunbiased_knowns = var(sqderr_linunbiased_knowns,0,2);
varse_linnjs = var(sqderr_linnjs,0,2);
varse_linnj_knowns = var(sqderr_linnj_knowns,0,2);
varse_ems = var(sqderr_ems,0,2);
varse_em_knowns = var(sqderr_em_knowns,0,2);
varse_blsslices = var(sqderr_blsslices,0,2);
if exist('tinv','file')
    c = tinv(0.95,ntests-1); % requires Statistics toolbox to compute
else
    c = 1.646380345428564; % for 95%, ntests = 1000
end
cilow_linunbiaseds = mse_linunbiaseds - (c/sqrt(ntests)).*sqrt(varse_linunbiaseds);
cihigh_linunbiaseds = mse_linunbiaseds + (c/sqrt(ntests)).*sqrt(varse_linunbiaseds);
cilow_linunbiased_knowns = mse_linunbiased_knowns - (c/sqrt(ntests)).*sqrt(varse_linunbiased_knowns);
cihigh_linunbiased_knowns = mse_linunbiased_knowns + (c/sqrt(ntests)).*sqrt(varse_linunbiased_knowns);
cilow_linnjs = mse_linnjs - (c/sqrt(ntests)).*sqrt(varse_linnjs);
cihigh_linnjs = mse_linnjs + (c/sqrt(ntests)).*sqrt(varse_linnjs);
cilow_linnj_knowns = mse_linnj_knowns - (c/sqrt(ntests)).*sqrt(varse_linnj_knowns);
cihigh_linnj_knowns = mse_linnj_knowns + (c/sqrt(ntests)).*sqrt(varse_linnj_knowns);
cilow_ems = mse_ems - (c/sqrt(ntests)).*sqrt(varse_ems);
cihigh_ems = mse_ems + (c/sqrt(ntests)).*sqrt(varse_ems);
cilow_em_knowns = mse_em_knowns - (c/sqrt(ntests)).*sqrt(varse_em_knowns);
cihigh_em_knowns = mse_em_knowns + (c/sqrt(ntests)).*sqrt(varse_em_knowns);
cilow_blsslices = mse_blsslices - (c/sqrt(ntests)).*sqrt(varse_blsslices);
cihigh_blsslices = mse_blsslices + (c/sqrt(ntests)).*sqrt(varse_blsslices);

clear *_temp;

% save results
if ~isempty(savepath)
    save(savepath);
end
% exit;
