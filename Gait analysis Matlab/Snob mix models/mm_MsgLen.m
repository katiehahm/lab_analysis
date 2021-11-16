function mm = mm_MsgLen(mm, data)

% Model and data structure
% ------------------------
n = size(data, 1);                  % Data set size
K = mm.nClasses;                    % Number of mixtures
a = mm.a;                           % Mixing proportions

% Assertion length for K (uniform prior)
% --------------------------------------
Ak = K*log(2);    % p(K) = 2^(-K)

% Assertion length for the mixing proportions
% -------------------------------------------
Aa = (K-1)*log(n)/2 - sum(log(a))/2 - gammaln(K);

% Detail length and assertion length for kn
% -----------------------------------------
% p = exp(-mm_Likelihood(mm, data, 1:mm.nModelTypes));
% p(p==0) = realmin; 
% r = bsxfun(@rdivide, p, sum(p,2));
p = -mm_Likelihood(mm, data, 1:mm.nModelTypes);
p = max(min(p, 700), -700); % log(p) \in [-700, 700]
ep = exp(p);
if(mm.nClasses > 1)
    r = bsxfun(@rdivide, ep, exp(logsumexp(p)));
else
    r = ones(mm.N, 1);
end
Nk = sum(r,1);      % how many things in each class

%An_L = -sum(log(sum(ep,2)));
An_L = -sum(logsumexp(p));
mm.L = An_L;

% Assertion length for the hyperparameters
% ----------------------------------------
Atheta = 0;
for i = 1:mm.nModelTypes
    switch mm.ModelTypes{i}.type
        %% Gaussian hyperparameters mu \in [mu0,mu1], tau \in [exp(-a),exp(+a)]
        case 'Gaussian'
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(2);
            end
            a_tau = FindPriorRange(tau);
            % two mu hyperparameters; each coded as log(n)/2
            % one tau hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);   
            
        %% Laplace hyperparameters mu \in [mu0,mu1], tau \in [exp(-a),exp(+a)]
        case 'Laplace'
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(2);
            end
            a_tau = FindPriorRange(tau);
            % two mu hyperparameters; each coded as log(n)/2
            % one tau hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);             
            
        %% Multivariate Gaussian
        case 'mvg'
            
            d = mm.ModelTypes{i}.nDim;
            Atheta = Atheta + sum(log(Nk)) * d;

        %% Single factor analysis model
        case 'sfa'
                       
            d = mm.ModelTypes{i}.nDim;
            sigma = zeros(K,d);
            for k =1:K
                sigma(k,:) = mm.class{k}.model{i}.theta(d+1:2*d);
            end
            a_sigma = FindPriorRange(sigma(:));            
            
            % two mu hyperparameters for each column of data; each coded as log(n)/2
            % one sigma hyperparameter coded as logstar(a_star)            
            Atheta = Atheta + sum(log(Nk))*d + K*d*log(2*a_sigma) + logstar(a_sigma);   
            
            
        %% Inverse Gaussian hyperparameters lambda \in [exp(-a), exp(+a)]
        case 'invGaussian'
            lambda = zeros(K,1);
            for k =1:K
                lambda(k) = mm.class{k}.model{i}.theta(2);
            end
            a_lambda = FindPriorRange(lambda);
            % mu hyperparameter; coded as log(n)/2
            % lambda hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk))/2 + K*log(2*a_lambda) + logstar(a_lambda);   
            
        %% Gaussian linear regression
        case 'linreg'     
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(1);
            end
            a_tau = FindPriorRange(tau);
            % two mu hyperparameters; each coded as log(n)/2
            % K region parameters coded as 1/2*log(n) each
            % tau hyperparameter coded as logstar(a)
            Atheta = Atheta + 3/2*sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);               
            
        %% Logistic regression
        case 'logreg'            
            Atheta = Atheta + 3/2*sum(log(Nk));
    end
end

% Assertion length for theta
% --------------------------
D = (K-1);      % number of mixing proportions
totalParams = (K-1);
nParams = 0;
for k = 1:K
    for i = 1:mm.nModelTypes

        model = mm.class{k}.model{i};           % model                
        switch model.type

            %% beta distribution
            case 'beta'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                ap = model.theta(1);
                bp = model.theta(2);  
                
                pgterm = psi(1,ap)*psi(1,bp) - (psi(1,ap)+psi(1,bp))*psi(1,ap+bp);
                %h_theta = -2*log(2) + 2*log(pi) + log1p(ap^2) + log1p(bp^2);
                % Let mu = a / (a+b), v = a + b
                % The prior below is obtained by assuming
                % mu ~ U(0,1), v ~ HC(2)
                % and transforming to (a,b) space
                h_theta = -log(4) + log(ap + bp) + log(pi) + log(4 + (ap + bp)^2);    
                F_theta = log(Nk(k)) + 0.5*log(pgterm);             
                
                AssLen = h_theta + F_theta;                  
            
            %% von Mises-Fisher distribution
            case 'vmf'
                nParams = mm.ModelTypes{i}.nDim;
                totalParams = totalParams + nParams;
                
                kappa = model.theta(1);
                %mu = model.theta(2:end);
                
                %% Prior densities
                d = nParams;
                if(d == 3)
                    logA = log(coth(kappa) - 1/kappa);
                elseif(d == 5)
                    logA = log(kappa/(-1+kappa*coth(kappa)) -(3/kappa));
                elseif(d == 7)
                    logA = log(kappa^3/(9+3*kappa^2-9*kappa*coth(kappa)) - 5/kappa - kappa/3);
                else
                    logA = log(besseli(d/2, kappa)) - log(besseli(d/2-1,kappa));
                end                
                
                h_kappa = -(d-1)*log(kappa) + (d+1)/2*log(1+kappa*kappa) - log(2) - gammaln((d+1)/2) + log(sqrt(pi)) + gammaln(d/2);
                h_mu = log(2) + (d/2)*log(pi) - gammaln(d/2);
                h_theta = h_mu + h_kappa;

                %% Fisher information
                logAp = log1p(- exp(2*logA) - (d-1)/kappa*exp(logA));
                F_theta = (d-1)*(log(Nk(k)) + log(kappa) + logA) + log(Nk(k)) + logAp;
                F_theta = F_theta * 0.5;                
                
                AssLen = h_theta + F_theta;  
            
            %% Weibull distribution
            case 'weibull'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                lambda = model.theta(1);
                k_wbl = model.theta(2);
                h_theta = -log(2) + log(pi) + log1p(lambda^2) -log(2) + log(pi) + log1p(k_wbl^2);
                F_theta = log(pi) - log(6)/2 - log(lambda) + log(Nk(k));
                AssLen = h_theta + F_theta;   
                
            %% Weibull distribution with type I fixed censoring
            case 'cfixweibull'                
                nParams = 2;
                totalParams = totalParams + nParams;
                                
                lambda = model.theta(1);
                k_wbl = model.theta(2);
                c = mm.ModelTypes{i}.c;    
                
                h_theta = -log(2) + log(pi) + log1p(lambda^2) -log(2) + log(pi) + log1p(k_wbl^2);
                p = 1 - exp(-(c/lambda)^k_wbl);
                logDetF = ((-0.144455E2)+(-0.194493E4).*p+(-0.188307E5).*p.^2+0.387201E5.* ...
                    p.^3+0.347085E5.*p.^4+(-0.115445E6).*p.^5+0.788731E5.*p.^6+( ...
                    -0.16047E5).*p.^7).*(1+0.237745E3.*p+0.48356E4.*p.^2+0.727136E4.* ...
                    p.^3+(-0.275255E5).*p.^4+0.11268E5.*p.^5+0.784586E4.*p.^6+( ...
                    -0.389537E4).*p.^7).^(-1);
                F_theta = log(Nk(k)) + 0.5*logDetF - log(lambda);
                AssLen = h_theta + F_theta;                   
                
            %% Laplace distribution
            case 'Laplace'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                b_lap = model.theta(2);   
                               
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                h_theta = log(Rmu) + log(b_lap);
                F_theta = log(Nk(k)) - 2*log(b_lap);
                AssLen = h_theta + F_theta;                   
            
            %% Univariate exponential
            case 'exp'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                lambda = model.theta;
                h_theta = -log(2) + log(pi) + log1p(lambda^2);
                F_theta = log(Nk(k))/2 - log(lambda);
                AssLen = h_theta + F_theta;
                
            %% Univariate exponential with type I random censoring
            case 'crndexp'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                alpha = model.theta(1); 
                beta = model.theta(2); 
                h_theta = -log(alpha) - log(beta) + log(alpha+beta) + 2*log(alpha*(beta+1)+beta);
                F_theta = log(Nk(k)) - 0.5*(log(beta) + log(alpha)) - log(alpha+beta);                
                AssLen = h_theta + F_theta;         
                
            %% Univariate exponential with fixed type I censoring                
            case 'cfixexp'
                nParams = 1;
                totalParams = totalParams + nParams;  
                
                c = mm.ModelTypes{i}.c;                
                theta = model.theta(1); 
                h_theta = -log(2) + log(pi) + log1p(theta*theta);
                F_theta = 0.5*log(Nk(k)) - log(theta) + 0.5*log(1-exp(-c/theta));    
                AssLen = h_theta + F_theta;                  
                
            %% Univariate gamma
            case 'gamma'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                mu = model.theta(1);
                phi = model.theta(2);
                
                h_theta = -log(2) + 2*log(pi) + log1p(mu^2) + log(phi)/2 + log1p(phi);
                F_theta =  log(Nk(k))- log(mu) + 1/2*log(phi*psi(1,phi) - 1);
                AssLen = h_theta + F_theta;
            
            %% Univariate k-nomial model
            case 'multi'
                
                % Hyperparameters
                alpha = mm.ModelTypes{i}.alpha;
                M = mm.ModelTypes{i}.nStates;
                A = mm.ModelTypes{i}.A;
                
                nParams = M-1;
                totalParams = totalParams + nParams;                
                theta = model.theta;
                
                h_theta = -gammaln(A) + sum(gammaln(alpha)) - sum( (alpha-1) .* log(theta) );
                F_theta = (M-1)/2*log(Nk(k)) - sum(log(theta))/2;
                
                AssLen = h_theta + F_theta;
            
            %% Single factor analysis model
            case 'sfa'
                % Parameters
                d = mm.ModelTypes{i}.nDim;                
                theta = model.theta;
                sigma = theta(d+1:2*d);
                if(~model.collapsed)    % correlated normal
                    nParams = 3*d + Nk(k);                    
                    a_sfa = theta(2*d+1:end);
                    beta = a_sfa ./ sigma;
                    b2 = beta'*beta;
                    v = model.v;
                               
                    Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;

                    Bk = zeros(d,1);
                    Bk(1:2)=[pi/2, pi];
                    for T=3:d
                        Bk(T)=2*pi * Bk(T-2)/(d-1);
                    end
                    Bk=Bk(d); 

                    F_theta = d/2*log(2*Nk(k)) + d/2*abs(log(Nk(k)*sum(v.^2)) - sum(v)^2) + (Nk(k)-2)/2*log(1+b2) - 3*sum(log(sigma));
                    h_theta = sum(log(Rmu)) + sum(log(sigma)) + (d+1)/2*log(1+b2) + log(Bk) + (Nk(k)/2)*log(2*pi) + v'*v/2;
                
                else    % uncorrelated model
                    nParams = 2*d;                    
                    
                    Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;

                    h_theta = sum(log(Rmu)) + sum(log(sigma));
                    F_theta = d*log(Nk(k)) - 2*sum(log(sigma)) + d*log(2)/2;                                       
                end
                
                totalParams = totalParams + nParams;                
                AssLen = h_theta + F_theta;                
                
            %% Multivariate Gaussian model
            case 'mvg'                
                d = mm.ModelTypes{i}.nDim;
                nParams = d + d*(d+1)/2;
                totalParams = totalParams + nParams;
                
                theta = model.theta;
                Sigma = reshape(theta(d+1:end),d,d);                
                
                R = cholcov(Sigma,0);              
                logDetSigma = 2*sum(log(diag(R)));
                
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                h_theta = sum(log(Rmu)) + (d+1)*logDetSigma + trace(R\(R'\eye(d)))/2 + log(2)*d*(d+1)/2 + logmvgamma(d,(d+1)/2);
                F_theta = d*(d+3)/4*log(Nk(k)) - d/2*log(2) -(d+2)/2*logDetSigma;
                AssLen = h_theta + F_theta;                
                                
            %% Univariate Gaussian model
            case 'Gaussian'
                nParams = 2;
                totalParams = totalParams + nParams;
                        
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                tau = model.theta(2);
                h_theta = log(Rmu) + log(tau);
                F_theta = log(Nk(k)) - 3*log(tau)/2 - log(2)/2;
                AssLen = h_theta + F_theta;
                
            %% Poisson model
            case 'Poisson'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                lambda = model.theta;
                h_theta = log(pi) + log(lambda)/2 + log1p(lambda);
                F_theta = log(Nk(k))/2 - log(lambda)/2;
                AssLen = h_theta + F_theta;
                
            %% Negative binomial
            case 'negb'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                mu = model.theta(1);
                phi = model.theta(2);
                
                h_theta = -2*log(2) + 2*log(pi) + log1p(phi^2) + log1p(mu^2);
                numerator = mu + phi*(mu + phi) * psi(1,phi)*((phi/(mu+phi))^phi - 1);
                F_theta = log(Nk(k)) - log(mu)/2 - log(mu+phi) + log(-numerator)/2;                
                AssLen = h_theta + F_theta;                
                
            %% geometric model
            case 'geometric'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                theta = model.theta;
                h_theta = -log(2-theta) + log(pi) + log(1-theta)/2 + log(theta^2 - theta + 1);
                F_theta = log(Nk(k))/2 - log(theta) - log(1-theta)/2;
                AssLen = h_theta + F_theta;
                
            %% Inverse Gaussian
            case 'invGaussian'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                % Parameters
                mu = model.theta(1);
                lambda = model.theta(2);                         
                
                % Hyperparameters
                mu0 = mm.ModelTypes{i}.mu0;        
                
                h_theta = -log(mu0)/2 + log(2) + 3*log(mu)/2 + log(lambda);
                F_theta = log(Nk(k)) - log(2)/2 - 3*log(mu)/2 - 3*log(lambda)/2;                
                AssLen = h_theta + F_theta;
                
                
            %% Gaussian linear regression
            case 'linreg'
                nParams = 2; % beta0 and tau; the others are handled within                
                CovIx = mm.ModelTypes{i}.CovIx; 
                Ivar  = mm.class{k}.model{i}.Ivar;
                totalParams = totalParams + nParams + length(CovIx);                
                
                
                P   = length(CovIx);
                ix  = ~isnan(data(:, Ivar)) & ~any(isnan(data(:,CovIx)),2);         
                X   = data(ix,mm.ModelTypes{i}.CovIx);                
                
                beta = model.theta(3:end);               
                logtau = log(model.theta(1));
                
                Xr = bsxfun(@times, X, sqrt(r(ix,k)));
                Kconst = beta'*(Xr'*Xr)*beta;
    
                log_kappa_P = -P*log(2) + log(P) + (1-P)*log(pi) + 2*psi(1)-P;
                logKmult = ( log_kappa_P + P * log(pi * Kconst) - 2*gammaln(P/2 + 1) );
                
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;               
                AssLen = log1p(exp(logKmult - P*logtau))/2 + log(Rmu) + logtau + log(Nk(k)) - log(2)/2 - 3*logtau/2;
                
            %% Logistic regression
            case 'logreg'
                nParams = 0; % all parameters (b0,b) coded within
                CovIx = mm.ModelTypes{i}.CovIx; 
                Ivar  = mm.class{k}.model{i}.Ivar;
                totalParams = totalParams + nParams + length(CovIx);        
                
                P   = length(CovIx);
                ix  = ~isnan(data(:, Ivar)) & ~any(isnan(data(:,CovIx)),2);         
                X   = data(ix,mm.ModelTypes{i}.CovIx);                
                Xr = bsxfun(@times, X, sqrt(r(ix,k)));
                
                beta0 = model.theta(1);
                beta  = model.theta(2:end);      
                
                % get mu
                lowerBnd = log(eps); 
                upperBnd = -lowerBnd;
                muLims = [eps, 1-eps];

                %% negative log-likelihood
                yhat = constrain(beta0 + X*beta, lowerBnd, upperBnd);
                mu = 1./(1 + exp(-yhat));

                if any(mu < muLims(1) | muLims(2) < mu)
                    mu = max(min(mu,muLims(2)),muLims(1));
                end                
                
                % assertion
                v = mu.*(1-mu);
                Xv = bsxfun(@times,Xr,sqrt(v));
                J = Xv'*Xv;
                yhat = Xr*beta;
                Kb = yhat'*yhat;

                logh = gammaln(P/2+1) + logdet(X'*X)/2 - (P/2)*log(pi) - (P/2)*log(Kb);
                log_kappa_P = -(P+1)*log(2) + log(P+1) + (1-(P+1))*log(pi) + 2*psi(1)-(P+1);
                logratio = log_kappa_P + logdet(J) - 2*logh;
                
                AssLen = log1p(exp(logratio))/2;                    

        end

        % Total assertion length and number of parameters
        Atheta = Atheta + AssLen;
        D = D + nParams;
    end
end

% Constant terms
% --------------
constant = mml_const(D) - gammaln(K+1);

% Total codelength for the mixture model 
% [see pp. 294, Statistical and Inductive Inference by Minimum Message Length]
% --------------------------------
mm.nParams = totalParams;   % total number of parameters
mm.Ak = Ak;                 % assertion length for K
mm.Aa = Aa;                 % assertion length for the class proportions
mm.Atheta = Atheta;         % assertion length for model parameters
mm.constant = constant;     % quantization constant
mm.msglen = Ak + Aa + Atheta + An_L + constant + log(n)/2;  % total codelength

% Compute AIC and BIC
% -------------------
mm.AIC = mm.L + (mm.nParams);
mm.BIC = mm.L + (mm.nParams/2)*log(n);

end