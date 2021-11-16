%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function mm = mm_EstimateTheta(mm, data, wClasses)
   
%% Estimate parameters
% For each class
for k = wClasses

    r = mm.r(:,k);  % posterior probabiity of item belonging to the class    
    
    % For each model (across columns), estimate parameters
    for i = 1:mm.nModelTypes
        
        model = mm.class{k}.model{i};           % model                
        y = data(:, mm.class{k}.model{i}.Ivar); % data 
        ix = ~isnan(y);                         % ignore missing data
        
        switch model.type
            
            %% Beta distribution
            case 'beta'
            
            Nk = sum(r(ix)); 
            S1 = sum(r(ix) .* log(y(ix)));     % sum log(y_i)
            S2 = sum(r(ix) .* log(1 - y(ix)));

            theta = fminunc(@(X) beta_msglen(Nk, S1, S2, X(1), X(2)), [0, 0], mm.opts.SearchOptions);
            model.theta = exp(theta);
            
            %% von Mises-Fisher distribution
            case 'vmf'
            
            ix  = ~any(isnan(y),2);
            kappa0 = 1;
            theta = fminunc(@(X) vmf_msglen(y(ix,:), r(ix), X(1)), kappa0, mm.opts.SearchOptions);
            
            S = sum(r(ix) .* y(ix,:));           % sum y_i
            S = S(:);
            model.theta = [exp(theta(1)); S / norm(S)];            
            
            %% Univariate Weibull distribution
            case 'weibull'
            theta = fminunc(@(X) weibull_msglen(y(ix), r(ix), X(1), X(2)), [log(mean(y(ix))), 0], mm.opts.SearchOptions);
            model.theta = max(1e-2, exp(theta(:)));
            
            %% Univariate Weibull distribution with fixed type I censoring
            case 'cfixweibull'
            c = mm.ModelTypes{i}.c;
            ix  = ~any(isnan(y),2);          
            
            k0 = 1;
            x0 = log( [ sum((1/(0.5 + sum(y(ix,2))))*y(ix,1).^k0).^(1/k0), k0] );            
            theta = fminunc(@(X) cfixweibull_msglen(y(ix,1), y(ix,2), c, r(ix), X(1), X(2)), x0, mm.opts.SearchOptions);
            model.theta = max(1e-2, exp(theta(:)));                
                            
            %% Univariate exponential
            case 'exp'
                
            % Sufficient statistics
            S  = sum(r(ix) .* y(ix));           % sum y_i

            % Estimate parameters
            Nk  = sum(r(ix));   
            Coeff = [-(Nk+1), S, 1-Nk, S];
            Coeff = Coeff ./ Coeff(1);            
            
            lambda = cubicroots(Coeff(2), Coeff(3), Coeff(4));
            
            model.theta = lambda;
            
            %% Exponential with random type I censoring
            case 'crndexp'
               
            ix  = ~any(isnan(y),2);                   
            
            % Sufficient statistics
            S = sum(r(ix) .* y(ix,1));           % sum y_i
            D = sum(r(ix) .* y(ix,2));           % sum delta_i
            
            Nk = sum(r(ix));
            
            % Estimates
            A = sqrt(2*(Nk+3)*S+(Nk-1)^2+S^2) - Nk + S + 1;
            alpha = A / (2*(Nk - D + 0.5));
            beta = A  / (2*(D + 0.5));
            model.theta = [ alpha; beta ];
            
            %% Exponential with fixed type I censoring
            case 'cfixexp'   
                
            ix  = ~any(isnan(y),2);                                        
            
            c = mm.ModelTypes{i}.c;
            S = sum(r(ix) .* y(ix,1));
            K = sum(r(ix) .* y(ix,2));
            
            x0 = log(S / (K + 0.5));
            theta = fminunc(@(X) cfixexp_msglen(S, K, c, X(1)), x0, mm.opts.SearchOptions);
            model.theta = exp(theta(:));                           
            
            %% Negative binomial distribution
            case 'negb'
            
            theta = fminunc(@(X) negb_msglen(y(ix), r(ix), X(1), X(2)), [0, 0], mm.opts.SearchOptions);
            model.theta = exp(theta(:));            
            
            %% Laplace distribution
            case 'Laplace'
            
            Nk = sum(r(ix));
            [mu,fval] = fminunc(@(MU) sum(r(ix) .* abs(y(ix) - MU)), 0, mm.opts.SearchOptions);
            %mu = medianw(y(ix),r(ix));
            %fval = sum(r(ix) .* abs(y(ix) - mu));
            b  = max(fval / (Nk - 1) , 1e-3 );
            
            model.theta = [mu; b];
            
            %% Gamma distribution
            case 'gamma'
                
            % Sufficient statistics
            S  = sum(r(ix) .* y(ix));
            L  = sum(r(ix) .* log(y(ix)));
            Nk = sum(r(ix)); 
            
%             % solve for phi
%             s_phi = log(S/Nk) - L/Nk ;
%             logphi_init = log(3 - s_phi + sqrt((s_phi - 3)^2 + 24*s_phi)) - log(12) - log(s_phi);
%             phi = exp( fminunc(@(X) gamma_msglen(Nk,S,L,X), logphi_init, mm.opts.SearchOptions) );
%             
%             % solve for mu
%             Coeff = [(1 + Nk*phi), -S*phi, (-1 + Nk*phi), -S*phi];
%             Coeff = Coeff ./ Coeff(1);                       
%             [x1,x2,x3] = cubicroots(Coeff(2), Coeff(3), Coeff(4));            
%             v = [x1,x2,x3];
%             mu=v(v>0);
%                 
%             model.theta = [mu; phi];

            [phi, mu, ~] = mm_EstimateGamma(S, L, Nk);
            model.theta = [mu; phi];
            
            %% Univariate k-nomial model
            case 'multi'
                
            % Hyperparameters
            alpha = mm.ModelTypes{i}.alpha;
            M = mm.ModelTypes{i}.nStates;
            A = mm.ModelTypes{i}.A;
            
            % Sufficient statistics
            Cnts = zeros(M, 1);
            R = r(ix);
            for j = 1:M
                Cnts(j) = sum(R(y(ix) == j));
            end
            Nk = sum(Cnts);

            % Estimate parameters
            model.theta = (Cnts + alpha - 1/2) ./ (Nk + A - M/2);
            
            %% Univariate Gaussian model
            case 'Gaussian'
        
            % Sufficient statistics
            s  = sum(r(ix) .* y(ix));           % sum y_i
            s2 = sum(r(ix) .* (y(ix).^2));      % sum y_i^2

            % Estimate parameters
            Nk  = sum(r(ix));                   % n
            mu  = s/Nk;
            tau = max( (s2 - 2*mu*s + Nk*mu^2) / (Nk-1), 1e-3);

            model.theta = [mu; tau];
            
            %% Single factor analysis
            case 'sfa'

            % Parameters
            ix  = ~any(isnan(y),2);        
            Nk  = sum(r(ix));                   % n            
            mu  = sum(bsxfun(@times, y(ix,:), r(ix))) / Nk;    
            xmu = bsxfun(@minus, y(ix,:), mu);
            rxmu = bsxfun(@times, xmu, sqrt(r(ix)));
            [a_sfa, sigma, v, collapsed] = mmlsfa(rxmu, Nk);
            
            model.theta = [mu(:); sigma(:); a_sfa(:)];   
            model.v     = v;
            model.collapsed = collapsed;

            
            %% Multivariate Gaussian model
            case 'mvg'

            % Estimate parameters
            d = mm.ModelTypes{i}.nDim;
            ix  = ~any(isnan(y),2);                                 
            Nk  = sum(r(ix));                   % n
            mu  = sum(bsxfun(@times, y(ix,:), r(ix))) / Nk;
            xmu = bsxfun(@minus, y(ix,:), mu);
            Sigma = (bsxfun(@times, xmu, r(ix))'*xmu + eye(d)) / (Nk+d);
            
            model.theta = [mu(:); Sigma(:)];                
                            
            %% Poisson model
            case 'Poisson'
                
            % Sufficient statistics
            s  = sum(r(ix) .* y(ix));           % sum y_i
            
            % Estimate parameters
            Nk  = sum(r(ix));                   % n            
            lambda = (sqrt(s^2 + 2*s*(Nk-1) + (Nk+1)^2) + s - Nk - 1) / 2 / Nk;
            
            model.theta = lambda;
            
            %% Geometric distribution
            case 'geometric'
                
            s  = sum(r(ix) .* y(ix));           % sum y_i
            Nk  = sum(r(ix));                   % n            
            geoR = roots( [Nk+s, 1-4*Nk-3*s, 1+6*Nk+3*s, -4-5*Nk-2*s, 2+2*Nk] );
            
            model.theta = geoR(imag(geoR) == 0 & geoR > 0 & geoR < 1);
            
            %% Inverse Gaussian
            case 'invGaussian'
                
            % Sufficient statstics
            S1 = sum(r(ix) .* y(ix));
            S2 = sum(r(ix) ./ y(ix));
            
            % Estimate parameters
            Nk = sum(r(ix));
            mu = S1 / Nk;
            lambda = max(1e-5, (S1*S2 - Nk^2) / (Nk - 1) / S1);
            
            model.theta = [mu; lambda];
            
            %% Linear regression
            case 'linreg'
                
            CovIx = mm.ModelTypes{i}.CovIx; 
            P   = length(CovIx);
            ix  = ix & ~any(isnan(data(:,CovIx)),2);         
            X   = data(ix, CovIx);
            y   = data(ix, mm.class{k}.model{i}.Ivar);

            % Estimate regression coefficients
            yr = y .* sqrt(r(ix));
            Xr = bsxfun(@times, X, sqrt(r(ix)));
            b = wridge(X, y, r(ix), 1);
            %b = lscov([ones(length(y),1), X], y, r(ix));

            % Estimate tau
            Nk  = sum(r(ix));   
            mu_z = Xr * b(2:end);
            Kconst = mu_z'*mu_z;
            e2 = sum( (yr - mu_z - sqrt(r(ix))*b(1)).^2 );
%            e2 = sum( (yr - mu_z - b(1)).^2 );
            log_kappa_P = -P*log(2) + log(P) + (1-P)*log(pi) + 2*psi(1)-P;
            logKmult = ( log_kappa_P + P * log(pi * Kconst) - 2*gammaln(P/2 + 1) );    
            
            logtau_init = log( max(1e-5, e2 / (Nk - P - 1)) );
            % uncomment below if actual MML estimate is required
            %logtau = fminunc(@(LOGTAU) linreg_msglentau(Nk, P, e2, logKmult, LOGTAU), logtau_init, opts.SearchOptions);
            logtau = logtau_init;
            tau = exp(max([logtau, -10]));
            
            model.theta = [tau; b];
            
            %% Logistic regression
            case 'logreg'
                
            CovIx = mm.ModelTypes{i}.CovIx; 
            P   = length(CovIx);
            ix  = ix & ~any(isnan(data(:,CovIx)),2);         
            X   = data(ix, CovIx);
            y   = data(ix, mm.class{k}.model{i}.Ivar);   
            
            % Estimate regression coefficients
            Xr = bsxfun(@times, X, sqrt(r(ix)));            
            yr = y .* sqrt(r(ix));
            
            b = fminunc(@(Z) logreg_msglen(r(ix), Xr, yr, X, y, Z(1), Z(2:end)), zeros(P+1,1), mm.opts.SearchOptions);            
            model.theta = b;
                        
        end
        
        mm.class{k}.model{i} = model;
        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = vmf_msglen(x, r, logk)

k = exp(logk);

%% Data
d = size(x,2);
n = sum(r);
R = sum(bsxfun(@times, x, r));
R = R(:);

%% Negative log-likelihood
if(d == 3)
    logbesseli = log(sinh(k)) + 0.5*log(2) - 0.5*logk - 0.5*log(pi);
    logA = log(coth(k) - 1/k);
else
    logbesseli = log(besseli(d/2-1,k));         % TODO numerical accuracy
    logA = log(besseli(d/2, k)) - logbesseli;
end

nll = -n*( (d/2-1)*log(k) - d/2*log(2*pi) - logbesseli ) - k*norm(R);

%% Prior densities
h_kappa = -(d-1)*logk + (d+1)/2*log(1+k*k) - log(2) - gammaln((d+1)/2) + log(sqrt(pi)) + gammaln(d/2);
h_mu = log(2) + (d/2)*log(pi) - gammaln(d/2);
h = h_mu + h_kappa;

%% Fisher information
logAp = log(1 - exp(2*logA) - (d-1)/k*exp(logA));
J = (d-1)*(log(n) + log(k) + logA) + log(n) + logAp;
J = J * 0.5;

%% Codelength
f = nll + h + J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = linreg_msglentau(n, p, e2, logKmult, logtau)

tau = exp(logtau);
f = n/2*logtau + e2/2/tau - logtau/2 + log1p(exp(logKmult - p*logtau))/2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = weibull_msglen(x, r, log_lambda, log_k)

k = exp(max(min(log_k, 1e2), -1e2));
lambda = exp(max(min(log_lambda, 1e2), -1e2));

F = -log_lambda;
h = log1p(lambda^2) + log1p(k^2);

%L = -log_k + k*log_lambda - (k-1)*log(x) + (x./lambda).^k;
z = x./lambda;
logz = log(z);
z2B = exp(k.*logz);
L = -((k-1).*logz + log(k./lambda)) + z2B;
%L(isinf(L)) = realmax;

f = F + h + r'*L;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msglen = cfixweibull_msglen(y, delta, c, r, log_lambda, log_k)

k = exp(min(max(log_k, -7), +7)); log_k = log(k);
lambda = exp(min(max(log_lambda, -7), +7)); log_lambda = log(lambda);

P = 2;
n = sum(r);

% assertion
kappa = 5/36/sqrt(3);
%h = -log(2) + log(pi) + log1p(lambda^2) + log(c) + 2*log(k/c) + c/k;
h = -log(2) + log(pi) + log1p(k^2) - log(2) + log(pi) + log1p(lambda^2);

%F = log(n) + log(det(wblconst_fisher(k, lambda, c)))/2;
p = 1 - exp(-(c./lambda).^k);
%logDetF = ((-0.107758E2)+(-0.878331E2).*p+0.194713E3.*p.^2+(-0.898213E2).* ...
%  p.^3).*(1+0.280711E2.*p+0.779886E1.*p.^2+(-0.240356E2).*p.^3).^(-1);
logDetF = ((-0.144455E2)+(-0.194493E4).*p+(-0.188307E5).*p.^2+0.387201E5.* ...
    p.^3+0.347085E5.*p.^4+(-0.115445E6).*p.^5+0.788731E5.*p.^6+( ...
    -0.16047E5).*p.^7).*(1+0.237745E3.*p+0.48356E4.*p.^2+0.727136E4.* ...
    p.^3+(-0.275255E5).*p.^4+0.11268E5.*p.^5+0.784586E4.*p.^6+( ...
    -0.389537E4).*p.^7).^(-1);
F = log(n) + 0.5*logDetF - log_lambda;
assertion = h + F + P/2*log(kappa);

% detail
z = y./lambda;
logz = log(z);
z2B = exp(k.*logz);
L = delta.*(-((k-1).*logz + log(k./lambda)) + z2B) + (1-delta).*(c/lambda)^k;    
L = r'*L;
detail = L + P/2;

msglen = assertion + detail;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = gamma_msglen(n, S, P, logphi)
% 
% phi = exp(logphi);
% Coeff = [(1 + n*phi), -S*phi, (-1 + n*phi), -S*phi];
% Coeff = Coeff ./ Coeff(1);                       
% [x1,x2,x3] = cubicroots(Coeff(2), Coeff(3), Coeff(4));
% v = [x1,x2,x3];
% mu=v(v>0);
% 
% f = n*gammaln(phi) + n*phi*log(mu/phi) - (phi-1)*P + phi/mu*S;  % neg ll
% f = f + (1/2)*log( phi*psi(1,phi) - 1 ) + log(1 + phi) + log(phi)/2;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, sigma, v, all_zeros] = mmlsfa(w, N)

%% setup
K = size(w,2);
V = w'*w;

% Initial guess for beta
[v,d] = eigs(V/N,1);
beta = v(:,end);
b2 = d(end);
beta = beta * sqrt(b2);

%% Estimate beta
done = false;
all_zeros = false;
reltol = 1e-4;
while(~done)
    b2 = beta'*beta;    
    
    if((N-1)*b2 < K)
       beta = zeros(K,1);
       sigma = std(w)';
       v = zeros(size(w,1),1);       
       all_zeros = true;
       done = true;
    else
       beta_old = beta;
       
       sigma = sqrt( sum(w.^2, 1)' / (N-1) ./ (1+beta.^2) );
       Y = V ./ (sigma * sigma');
       beta = Y*beta * (1 - K/(N-1)/b2) / (N-1) / (1+b2);      
              
       %done = norm(beta - beta_old) < 1e-3;
       if(norm((beta - beta_old) ./ (1.0 + abs(beta_old)), Inf) < reltol)
            done = true;
       end         
    end
end

%% Other parameters
b2 = beta'*beta;
a = sigma .* beta;
if(~all_zeros)
    y = bsxfun(@rdivide, w, sigma');    
    v = (y * beta) * (1 - K/(N-1)/b2) / (1 + b2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = cfixexp_msglen(S, K, c, log_theta)

theta = exp(log_theta);
L = K*log_theta + S/theta;
%L = r'*(delta.*log_theta + y/theta);
F = - log_theta + 0.5*log(1-exp(-c/theta));
h = log1p(theta*theta);

f = L + F + h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = negb_msglen(x, r, logmu, logphi)

n = sum(r);
mu = exp(logmu);
phi = exp(logphi);

L = -gammaln(x+phi) + gammaln(phi) + gammaln(x+1) - x .* log(mu/(mu+phi)) - phi*log(phi/(mu+phi));
h = -2*log(2) + 2*log(pi) + log1p(phi^2) + log1p(mu^2);

numerator = mu + phi*(mu + phi) * psi(1,phi)*((phi/(mu+phi))^phi - 1);
J = log(n) - log(mu)/2 - log(mu+phi) + log(-numerator)/2;

f = r'*L + h + J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = beta_msglen(Nk, S1, S2, loga, logb)

a = exp(loga);
b = exp(logb);

L = -(a-1)*S1 - (b-1)*S2 + Nk*betaln(a,b);

pgterm = psi(1,a)*psi(1,b) - (psi(1,a)+psi(1,b))*psi(1,a+b);
%h = -2*log(2) + 2*log(pi) + log1p(a^2) + log1p(b^2);
h = -log(4) + log(a + b) + log(pi) + log(4 + (a + b)^2);
J = log(Nk) + 0.5*log(pgterm);

f = h + J + L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msglen = logreg_msglen(r, Xr, yr, X, y, b0, b)

% number of parameters
P = length(b);
K = P + 1;

% negative log-likelihood
[L,mu] = logregnll(r,X,y,b0,b);

% assertion
v = mu.*(1-mu);
Xv = bsxfun(@times,Xr,sqrt(v));
J = Xv'*Xv;
yhat = Xr*b;
Kb = yhat'*yhat;
logh = gammaln(P/2+1) + logdet(Xr'*Xr)/2 - (P/2)*log(pi) - (P/2)*log(Kb);
log_kappa_K = -K*log(2) + log(K) + (1-K)*log(pi) + 2*psi(1) - K;
logratio = log_kappa_K + logdet(J) - 2*logh;

assertion = log1p(exp(logratio))/2;

% msglen
msglen = assertion + L + K/2;

end

function [f,mu,yhat] = logregnll(r,X,y,b0,b)

lowerBnd = log(eps); 
upperBnd = -lowerBnd;
muLims = [eps, 1-eps];

%% negative log-likelihood
yhat=constrain(b0 + X*b, lowerBnd, upperBnd);
mu=1./(1 + exp(-yhat));

if any(mu < muLims(1) | muLims(2) < mu)
    mu = max(min(mu,muLims(2)),muLims(1));
end

f = -(y.*log(mu) + (1.0-y).*log(1.0-mu));
f = r'*f;


end