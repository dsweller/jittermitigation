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

function varargout = gibbsslicesampler(K, M, y, h, sigma_x, sigma_w, sigma_z, varargin)

% initialization
N = length(y);
ns_minus_ks = crossdiff((0:N-1)'./M,(0:K-1));
z0 = zeros(N,1);
x0 = zeros(K,1);
if isempty(sigma_x), sigma_x = 1; end
if isempty(sigma_w), sigma_w = 0.1; end
if isempty(sigma_z), sigma_z = 0.1; end
sigma_x0 = reshape(sigma_x,1,[]);
sigma_w0 = reshape(sigma_w,1,[]);
sigma_z0 = reshape(sigma_z,1,[]);
sigma_xrandom = 0;
sigma_wrandom = 0;
sigma_zrandom = 0;
sigma_xparams = {};
sigma_wparams = {};
sigma_zparams = {};
shrinkagemethod = 'reject';
plotshrinkmethod = false;
chains = 1;

numiters = 75;
numbootstrap = 25;
maxshrinkiters = 200;
shrinkagethreshold = 25; % on log scale, only used for rejection-midpoint-threshold shrinkage algorithm

for ia = 1:2:nargin-8
    argkey = varargin{ia};
    argval = varargin{ia+1};
    switch lower(argkey)
        case 'chains'
            chains = argval;
        case 'z0'
            z0 = argval;
        case 'x0'
            x0 = argval;
        case 'randomsigmax'
            sigma_xrandom = 1;
            sigma_xparams = argval;
        case 'randomsigmaw'
            sigma_wrandom = 1;
            sigma_wparams = argval;
        case 'randomsigmaz'
            sigma_zrandom = 1;
            sigma_zparams = argval;
        case 'sigma_x0'
            sigma_x0 = argval;
        case 'sigma_w0'
            sigma_w0 = argval;
        case 'sigma_z0'
            sigma_z0 = argval;
        case 'numiters'
            numiters = argval;
        case 'numbootstrap'
            numbootstrap = argval;
        case 'shrinkagemethod'
            shrinkagemethod = argval;
        case 'maxshrinkiters'
            maxshrinkiters = argval;
        case 'shrinkagethreshold'
            shrinkagethreshold = argval;
        case 'plotshrinkmethod'
            plotshrinkmethod = argval;
        otherwise
            warning('STIR:invalidArg','Invalid argument \"%s\".',argkey);
    end
end

chains = [chains,size(z0,2),size(x0,2)];
if sigma_xrandom, chains = [chains,length(sigma_x0)]; end
if sigma_wrandom, chains = [chains,length(sigma_w0)]; end
if sigma_zrandom, chains = [chains,length(sigma_z0)]; end
numchains = max(chains);
if size(z0,2) < numchains, z0 = z0(:,1)*ones(1,numchains); end
if size(x0,2) < numchains, x0 = x0(:,1)*ones(1,numchains); end
if sigma_xrandom && size(sigma_x0,2) < numchains, sigma_x0 = sigma_x0(1)*ones(1,numchains); end
if sigma_wrandom && size(sigma_w0,2) < numchains, sigma_w0 = sigma_w0(1)*ones(1,numchains); end
if sigma_zrandom && size(sigma_z0,2) < numchains, sigma_z0 = sigma_z0(1)*ones(1,numchains); end

z_hat = z0;
x_hat = x0;

sigma_x_hat = zeros(1,numchains)+sigma_x0; % in non-random case, equals sigma_x
sigma_w_hat = zeros(1,numchains)+sigma_w0; % in non-random case, equals sigma_w
sigma_z_hat = zeros(1,numchains)+sigma_z0; % in non-random case, equals sigma_z

z_hats = zeros(N,numchains,numiters+1);
x_hats = zeros(K,numchains,numiters+1);
all_vars = [z_hat;x_hat];
if sigma_xrandom
    sigma_x_hats = zeros(1,numchains,numiters+1);
    all_vars = [all_vars;sigma_x_hat];
end
if sigma_wrandom
    sigma_w_hats = zeros(1,numchains,numiters+1);
    all_vars = [all_vars;sigma_w_hat];
end
if sigma_zrandom
    sigma_z_hats = zeros(1,numchains,numiters+1);
    all_vars = [all_vars;sigma_z_hat];
end
all_vars_sum = all_vars;
all_vars_sqd_sum = zeros([size(all_vars,1),size(all_vars,1),numchains]);
for m=1:numchains
    all_vars_sqd_sum(:,:,m) = all_vars_sqd_sum(:,:,m) + all_vars(:,m)*all_vars(:,m)';
end

cov_interchain = zeros([N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,numbootstrap+numiters]);
cov_intrachain = zeros([N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,numbootstrap+numiters]);
Vhat = zeros(N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,N+K+sigma_xrandom+sigma_wrandom+sigma_zrandom,numbootstrap+numiters);
PSRF = zeros(1,numbootstrap+numiters);

ns = repmat((0:N-1)',[1,numchains]);
nchains = repmat((1:numchains),N,1);
ns_minus_ks = repmat(ns_minus_ks,[1,1,numchains]);
shrinkiters = zeros(N,numchains,numiters+numbootstrap);

for ii = 1-numbootstrap:numiters
    
    if ii == 1 % set initial values
        z_hats(:,:,1) = z_hat;
        x_hats(:,:,1) = x_hat;
        if sigma_xrandom, sigma_x_hats(:,:,1) = sigma_x_hat; end
        if sigma_wrandom, sigma_w_hats(:,:,1) = sigma_w_hat; end
        if sigma_zrandom, sigma_z_hats(:,:,1) = sigma_z_hat; end
    end
    
    % perform Gibbs sampling iteration
    
    % generate z_hat using slice sampling

    % compute p(z_n|y_n,x,sigma_z,sigma_w)
    sigma_w_hat_rep = repmat(sigma_w_hat,N,1);
    sigma_z_hat_rep = repmat(sigma_z_hat,N,1);
    pz = pznmixnorm(z_hat(:), y(ns(:)+1), ns(:), h, x_hat(:,nchains), K, M, sigma_w_hat_rep(:), sigma_z_hat_rep(:));

    % compute u
    us = pz + log(rand(N*numchains,1));
    
    % assemble initial slice interval
    Rs = sigma_z_hat_rep(:).*sqrt(-2*us);
    Ls = -Rs;
    
    if plotshrinkmethod && ii == 1 && usejava('awt')
        hplotshrink = figure;
        subplot(2,1,1);
        ncenter = M*floor(K/2)+1;
        z_hat_center = z_hat(ncenter,1);
        z_hat_grid = linspace(Ls(ncenter,1),Rs(ncenter,1),200).';
        pz_grid = pznmixnorm(z_hat_grid,y(ncenter,1), ncenter-1, h, repmat(x_hat(:,1),[1,length(z_hat_grid)]), K, M, sigma_w_hat(1), sigma_z_hat(1));
        plot(z_hat_grid,log(exp(pz_grid)),'k'); hold on;
        % create filled region
        ptstart = [];
        for ipt = 1:length(z_hat_grid)
            if pz_grid(ipt) >= us(ncenter,1) && isempty(ptstart), ptstart = ipt; end
            if pz_grid(ipt) < us(ncenter,1) && ~isempty(ptstart)
                polyx = z_hat_grid([ptstart,ptstart:ipt-1,ipt-1]);
                polyy = [us(ncenter,1);pz_grid(ptstart:ipt-1);us(ncenter,1)];
                patch('XData',polyx,'YData',polyy,'FaceColor',0.9.*[1,1,1]);
                ptstart = [];
            end
        end
        if ~isempty(ptstart)
            polyx = z_hat_grid([ptstart,ptstart:length(z_hat_grid),length(z_hat_grid)]);
            polyy = [us(ncenter,1);pz_grid(ptstart:end);us(ncenter,1)];
            patch('XData',polyx,'YData',polyy,'FaceColor',0.9.*[1,1,1]);
        end
        line([Ls(ncenter,1);Rs(ncenter,1)],[us(ncenter,1);us(ncenter,1)],'LineStyle',':','Color','k');
        text(Rs(ncenter,1),us(ncenter,1),'slice','HorizontalAlignment','right','VerticalAlignment','bottom','Margin',2);
        xlim([Ls(ncenter,1),Rs(ncenter,1)]);
        if strcmpi(shrinkagemethod,'rej-midpt-th')
            line([Ls(ncenter,1);Rs(ncenter,1)],[us(ncenter,1);us(ncenter,1)]-shrinkagethreshold,'LineStyle',':','Color','k');
            text(Rs(ncenter,1),us(ncenter,1)-shrinkagethreshold,'midpt. threshold','HorizontalAlignment','right','VerticalAlignment','top','Margin',2);
        end
        ylabel('ptilde');
        ylims = ylim();
        if ylims(1) < us(ncenter,1)-100, ylims(1) = us(ncenter,1)-100; ylim(ylims); end
        line([z_hat_center;z_hat_center],[ylims(1);pz(ncenter,1)],'LineStyle',':','Color','k');
        set(gca,'Box','off');
        subplot(2,1,2);
        xlim([Ls(ncenter,1),Rs(ncenter,1)]);
        ylabel('# shrinkage iterations');
        set(gcf,'Color','w','ResizeFcn',@shrinkplot_resizefcn);
    else
        hplotshrink = [];
    end
    
    % perform shrinkage iterations
    notfound = true(N,numchains);
    nowshrinkiters = zeros(N,numchains);
    numnotfound = N*numchains;
    for is = 1:maxshrinkiters
        
        if plotshrinkmethod && ii == 1 && notfound(ncenter,1)
            nindex = 1:N; nindex = nindex(notfound(:,1)); nindex = find(nindex==ncenter,1);
        end
        
        % generate candidate
        zproposed = rand(numnotfound,1).*(Rs(notfound)-Ls(notfound))+Ls(notfound);
        pzproposed = pznmixnorm(zproposed, y(ns(notfound)+1), ns(notfound), h, x_hat(:,nchains(notfound)), K, M, sigma_w_hat_rep(notfound), sigma_z_hat_rep(notfound));
        
        % evaluate candidates
        found = (pzproposed - us(notfound)) > 0; % don't treat as an add op
        foundnow = false(N,numchains);
        foundnow(notfound) = found;
        z_hat(foundnow) = zproposed(found);
        nowshrinkiters(foundnow) = is;
        notfound(foundnow) = false;
        numnotfound = numnotfound - sum(found(:));
        
        if plotshrinkmethod && ii == 1 && ~isempty(hplotshrink) && ishandle(hplotshrink) && usejava('awt') % ignore for op count
            figure(hplotshrink);
            axs = get(hplotshrink,'Children');
            set(hplotshrink,'CurrentAxes',axs(1));
            line([Ls(ncenter,1),Ls(ncenter,1),Rs(ncenter,1);Rs(ncenter,1),Ls(ncenter,1),Rs(ncenter,1)],[is,is-0.125,is-0.125;is,is+0.125,is+0.125],'Color','k');
            xlims = xlim();
            text(xlims(1),is,sprintf('%d:',is),'HorizontalAlignment','right','Margin',4);
            if ~isempty(nindex)
                if found(nindex)
                    line([zproposed(nindex);zproposed(nindex)],[is;is],'LineStyle','none','Marker','o','MarkerEdgeColor','k');
                    text(zproposed(nindex),is,'z_n^{(i)} \leftarrow z','HorizontalAlignment','center','VerticalAlignment','top','Margin',2);
                else
                    line([zproposed(nindex);zproposed(nindex)],[is;is],'LineStyle','none','Marker','x','MarkerEdgeColor','k','MarkerSize',8);
                    if zproposed(nindex) < z_hat_center
                        text(zproposed(nindex),is,'L \leftarrow z','HorizontalAlignment','center','VerticalAlignment','top','Margin',2);
                    else
                        text(zproposed(nindex),is,'R \leftarrow z','HorizontalAlignment','center','VerticalAlignment','top','Margin',2);
                    end
                end
            end
            
            if ~notfound(ncenter,1)
                line([z_hat_center;z_hat_center],[0.5;is+0.5],'LineStyle',':','Color','k');
                ylim([0.5,is+0.5]);
                text(z_hat_center,1,'z_n^{(i-1)}','HorizontalAlignment','left','VerticalAlignment','top','Margin',2);
                pos = get(axs(2),'Position'); posti = get(axs(2),'TightInset');
                newpos = get(axs(1),'Position'); newpos(4) = pos(2)-posti(2)-newpos(2);
                set(axs(1),'Position',newpos,'ActivePositionProperty','position','YDir','reverse','Box','off','YTick',1:is,'Visible','off');
                set(get(axs(1),'YLabel'),'Visible','on');
                hplotshrink = [];
            end
        end
        
        if numnotfound == 0
            break;
        end

        belowthresholdfull = notfound;
        if strcmpi(shrinkagemethod,'rej-midpt-th')
            belowthresholdfull(notfound) = (pzproposed(~found) - us(notfound)) < -shrinkagethreshold;
        end
        % shrink intervals
        switch lower(shrinkagemethod)
            case {'reject','rej-midpt','rej-midpt-th'}
                shrinkLs = false(N,numchains);
                shrinkRs = false(N,numchains);
                zrejected = zproposed(~found);
                zLs = zrejected < z_hat(notfound);
                zRs = ~zLs;
                shrinkLs(notfound) = zLs;
                shrinkRs(notfound) = zRs;
                Ls(shrinkLs) = zrejected(zLs);
                Rs(shrinkRs) = zrejected(zRs);
        end
        switch lower(shrinkagemethod)
            case {'midpoint','rej-midpt','rej-midpt-th'}
                shrinkLs = false(N,numchains);
                shrinkRs = false(N,numchains);
                zrejected = (Ls(belowthresholdfull)+Rs(belowthresholdfull))./2;
                zLs = zrejected < z_hat(belowthresholdfull);
                zRs = ~zLs;
                shrinkLs(belowthresholdfull) = zLs;
                shrinkRs(belowthresholdfull) = zRs;
                Ls(shrinkLs) = zrejected(zLs);
                Rs(shrinkRs) = zrejected(zRs);
                if plotshrinkmethod && ii == 1 && belowthresholdfull(ncenter,1) && ~isempty(hplotshrink) && ishandle(hplotshrink) && usejava('awt')
                    if shrinkLs(ncenter,1)
                        text(Ls(ncenter,1),is+0.5,'L \leftarrow midpt.','HorizontalAlignment','center','VerticalAlignment','middle','Margin',2);
                    else
                        text(Rs(ncenter,1),is+0.5,'R \leftarrow midpt.','HorizontalAlignment','center','VerticalAlignment','middle','Margin',2);
                    end
                end
        end
    end
    shrinkiters(:,:,ii+numbootstrap) = nowshrinkiters;
    
    if ii >= 1, z_hats(:,:,ii+1) = z_hat; end
    all_vars(1:N,:) = z_hat;

    % generate x_hat from normal distribution
    ts = ns_minus_ks + permute(z_hat(:,:,ones(1,K)),[1,3,2]);
    Hz = h(ts,K);
    for ic = 1:numchains
        HTH = Hz(:,:,ic).'*Hz(:,:,ic);
        HTH = HTH + (sigma_w_hat(ic)./sigma_x_hat(ic)).^2.*eye(K);
        [U,S,V] = svd(HTH,'econ'); s = diag(S);
        
        mu_x = V*diag(s.^-1)*U'*(Hz(:,:,ic).'*y);
        x_hat(:,ic) = sigma_w_hat(ic).*V*diag(s.^(-0.5))*U'*randn([K,1])+mu_x;
    end

    if ii >= 1, x_hats(:,:,ii+1) = x_hat; end
    all_vars(N+(1:K),:) = x_hat;
    
    if sigma_xrandom
        % generate sigma_w_hat using inverse gamma
        alpha = zeros(1,numchains) + (sigma_xparams{1} + K/2);
        beta = zeros(1,numchains) + sigma_xparams{2};
        for ic = 1:numchains
            beta(ic) = beta(ic) + sum(x_hat(:,ic).^2)/2;
        end
        gamresult = gammar(alpha);
        sigma_x_hat = sqrt(beta./gamresult);
    
        if ii >= 1, sigma_x_hats(:,:,ii+1) = sigma_x_hat; end
        all_vars(N+K+1,:) = sigma_x_hat;
    end
    
    if sigma_wrandom
        % generate sigma_w_hat using inverse gamma
        alpha = zeros(1,numchains) + (sigma_wparams{1} + N/2);
        beta = zeros(1,numchains) + sigma_wparams{2};
        for ic = 1:numchains
            beta(ic) = beta(ic) + sum((y - Hz(:,:,ic)*x_hat(:,ic)).^2)/2;
        end
        gamresult = gammar(alpha);
        sigma_w_hat = sqrt(beta./gamresult);
    
        if ii >= 1, sigma_w_hats(:,:,ii+1) = sigma_w_hat; end
        all_vars(N+K+sigma_xrandom+1,:) = sigma_w_hat;
    end
    
    if sigma_zrandom
        % generate sigma_z_hat using inverse gamma
        alpha = zeros(1,numchains) + (sigma_zparams{1} + N/2);
        beta = zeros(1,numchains) + sigma_zparams{2};
        for ic = 1:numchains
            beta(ic) = beta(ic) + sum(z_hat(:,ic).^2)/2;
        end
        gamresult = gammar(alpha);
        sigma_z_hat = sqrt(beta./gamresult);
        
        if ii >= 1, sigma_z_hats(:,:,ii+1) = sigma_z_hat; end
        all_vars(N+K+sigma_xrandom+sigma_wrandom+1,:) = sigma_z_hat;
    end

    all_vars_sum = all_vars_sum + all_vars;
    for m=1:numchains
        all_vars_sqd_sum(:,:,m) = all_vars_sqd_sum(:,:,m) + all_vars(:,m)*all_vars(:,m)';
    end
    if numchains > 1 && ii >=1-numbootstrap && nargout >= 6+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom
        all_vars_sum_sqd = zeros([size(all_vars,1),size(all_vars,1),numchains]);
        for m=1:numchains
            all_vars_sum_sqd(:,:,m) = all_vars_sum(:,m)*all_vars_sum(:,m)';
        end
        % compute covariances
        cov_interchain(:,:,ii+numbootstrap) = cov(all_vars_sum.'./(ii+numbootstrap+1));
        cov_intrachain(:,:,ii+numbootstrap) = sum(all_vars_sqd_sum-all_vars_sum_sqd./(ii+numbootstrap+1),3)./(numchains*(ii+numbootstrap));
        
        % compute Vhat and PSRF
        Vhat(:,:,ii+numbootstrap) = ((ii+numbootstrap)/(ii+numbootstrap+1)).*cov_intrachain(:,:,ii+numbootstrap)+((numchains+1)/numchains).*cov_interchain(:,:,ii+numbootstrap);
        PSRF(ii+numbootstrap) = (ii+numbootstrap)/(ii+numbootstrap+1)+(numchains+1)/numchains*norm((cov_intrachain(:,:,ii+numbootstrap)+1e-12*eye(size(all_vars,1)))\cov_interchain(:,:,ii+numbootstrap),2);
    end
end

% compute estimates
z_est = mean(z_hats(:,:,2:end),3);
x_est = mean(x_hats(:,:,2:end),3);
if sigma_xrandom, sigma_x_est = mean(sigma_x_hats(:,:,2:end),3); end
if sigma_wrandom, sigma_w_est = mean(sigma_w_hats(:,:,2:end),3); end
if sigma_zrandom, sigma_z_est = mean(sigma_z_hats(:,:,2:end),3); end

if nargout >= 1, varargout{1} = z_est; end
if nargout >= 2, varargout{2} = x_est; end
if sigma_xrandom && nargout >= 3, varargout{3} = sigma_x_est; end
if sigma_wrandom && nargout >= 3+sigma_xrandom, varargout{3+sigma_xrandom} = sigma_w_est; end
if sigma_zrandom && nargout >= 3+sigma_xrandom+sigma_wrandom, varargout{3+sigma_xrandom+sigma_wrandom} = sigma_z_est; end
if nargout >= 3+sigma_xrandom+sigma_wrandom+sigma_zrandom, varargout{3+sigma_xrandom+sigma_wrandom+sigma_zrandom} = z_hats; end
if nargout >= 4+sigma_xrandom+sigma_wrandom+sigma_zrandom, varargout{4+sigma_xrandom+sigma_wrandom+sigma_zrandom} = x_hats; end
if sigma_xrandom && nargout >= 5+1+sigma_wrandom+sigma_zrandom, varargout{5+1+sigma_wrandom+sigma_zrandom} = sigma_x_hats; end
if sigma_wrandom && nargout >= 5+2*sigma_xrandom+1+sigma_zrandom, varargout{5+2*sigma_xrandom+1+sigma_zrandom} = sigma_w_hats; end
if sigma_zrandom && nargout >= 5+2*sigma_xrandom+2*sigma_wrandom+1, varargout{5+2*sigma_xrandom+2*sigma_wrandom+1} = sigma_z_hats; end
if nargout >= 5+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom, varargout{5+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom} = shrinkiters; end
if nargout >= 6+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom, varargout{6+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom} = cov_interchain; end
if nargout >= 7+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom, varargout{7+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom} = cov_intrachain; end
if nargout >= 8+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom, varargout{8+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom} = Vhat; end
if nargout >= 9+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom, varargout{9+2*sigma_xrandom+2*sigma_wrandom+2*sigma_zrandom} = PSRF; end

end

function pz = pznmixnorm(z, y, ns, h, x, K, M, sigma_w, sigma_z)

ks = (0:K-1);
h_n = h(crossdiff(ns./M+z,ks),K);
pz = -(y-sum(h_n.*x.',2)).^2./(2.*sigma_w.^2)-z.^2./(2.*sigma_z.^2);

end

function shrinkplot_resizefcn(src,evt) %#ok<INUSD>

axs = get(gcf,'Children');
pos = get(axs(2),'Position'); posti = get(axs(2),'TightInset');
newpos = get(axs(1),'Position'); newpos(4) = pos(2)-posti(2)-newpos(2);
set(axs(1),'Position',newpos);

end
