function L = mm_Likelihood(mm, Y, wModels)

n = size(Y,1);
K = mm.nClasses;   % number of mixtures
a = mm.a;          % mixing proportions

L = zeros(n,K);

%% Get likelihoods
% For each class
for k = 1:K
    % For each model (across columns), get likelihoods
    for i = wModels
        
        subL = zeros(n, 1);
                    
        m = mm.class{k}.model{i};            
        I = ~any(isnan( Y(:,m.Ivar) ), 2);
        
        switch m.type            
            
            %% beta distribution
            case 'beta'
                ap = m.theta(1);
                bp = m.theta(2);
                
                subL(I) = betaln(ap,bp) - (ap-1)*log(Y(I, m.Ivar)) - (bp-1)*log(1-Y(I, m.Ivar));
            
            %% von Mises-Fisher 
            case 'vmf' 
                kappa = m.theta(1);
                mu = m.theta(2:end);
                
                d = length(m.Ivar);
                if(d == 3)
                    logbesseli = log(sinh(kappa)) + 0.5*log(2) - 0.5*log(kappa) - 0.5*log(pi);
                else
                    % TODO: improve numerical accuracy
                    logbesseli = log(besseli(d/2-1,kappa));
                end
                
                subL(I) = -( (d/2-1)*log(kappa) - d/2*log(2*pi) - logbesseli ) - kappa*Y(I, m.Ivar)*mu;
            
            %% Weibull
            case 'weibull'
                lambda = m.theta(1);
                k_wbl = m.theta(2);                
                y = Y(I, m.Ivar(1)); 
                
                z = y./lambda;
                logz = log(z);
                z2B = exp(k_wbl.*logz);
                subL(I) = -((k_wbl-1).*logz + log(k_wbl./lambda)) + z2B;
                
            %% Weibull with fixed type I censoring
            case 'cfixweibull'
                lambda = m.theta(1);
                k_wbl = m.theta(2);
                c = mm.ModelTypes{i}.c;
                
                y = Y(I, m.Ivar(1));
                delta = Y(I, m.Ivar(2));    
                
                z = y./lambda;
                logz = log(z);
                z2B = exp(k_wbl.*logz);
                subL(I) = delta.*(-((k_wbl-1).*logz + log(k_wbl./lambda)) + z2B) + (1-delta).*(c/lambda)^k_wbl;                
            
            %% Negative binomial
            case 'negb'
                mu = m.theta(1);
                phi = m.theta(2);
                
                subL(I) = -gammaln(Y(I, m.Ivar)+phi) + gammaln(phi) + gammaln(Y(I, m.Ivar)+1) - Y(I, m.Ivar).*log(mu/(mu+phi)) - phi*log(phi/(mu+phi));
                
            %% exponential
            case 'exp'                
                lambda = m.theta;
                subL(I) = log(lambda) +  Y(I, m.Ivar) ./ lambda;
                
            %% exponential with type I random censoring
            case 'crndexp'
                alpha = m.theta(1);
                beta = m.theta(2);
                y = Y(I, m.Ivar(1));
                D = Y(I, m.Ivar(2));                
                subL(I) = D*log(beta) + (1-D)*log(alpha) + (1/alpha+1/beta)*y;
                
            %% exponential with type I fixed censoring
            case 'cfixexp'
                theta = m.theta(1);
                y = Y(I, m.Ivar(1));
                D = Y(I, m.Ivar(2));                
                subL(I) = D.*log(theta) + y./theta;
                
            %% Laplace
            case 'Laplace'
                mu = m.theta(1);
                b_lap = m.theta(2);
                
                subL(I) = log(2*b_lap) + abs(Y(I, m.Ivar) - mu)/b_lap;
                
            %% gamma
            case 'gamma'
                mu = m.theta(1);
                phi = m.theta(2);
                
                subL(I) = gammaln(phi) + phi*log(mu/phi) - (phi-1)*log(Y(I, m.Ivar)) + phi/mu*Y(I, m.Ivar);
            
            %% Univariate k-nomial
            case 'multi'
                
                % Parameters
                theta = m.theta;
                subL(I) = -log( theta(Y(I, m.Ivar)) );
                
            %% Single factor analysis model
            case 'sfa'
                
                % Parameters
                d = mm.ModelTypes{i}.nDim;
                theta = m.theta;
                mu = theta(1:d);
                sigma = theta(d+1:2*d);
                a_sfa = theta(2*d+1:end);
                
                x0 = bsxfun(@minus, Y(I, m.Ivar), mu');
                vrep = repmat(m.v, [1,d]);
                w = x0 - bsxfun(@times,vrep,a_sfa');
                e2 = sum(bsxfun(@rdivide, w.^2, sigma.^2'),2);
                subL(I) = (d/2)*log(2*pi) + sum(log(sigma)) + e2/2;
                
                
            %% Multivariate Gaussian model
            case 'mvg'
                
                % Parameters
                d = mm.ModelTypes{i}.nDim;
                theta = m.theta;
                mu = theta(1:d);
                Sigma = reshape(theta(d+1:end),d,d);
                
                X0 = bsxfun(@minus, Y(I, m.Ivar), mu');
                [R, err] = cholcov(Sigma,0);              
                logSqrtDetSigma = sum(log(diag(R)));
                xRinv = X0 / R;
                quadform = sum(xRinv.^2, 2);
                
                % Negative log-likelihood
                subL(I) = 0.5*quadform + logSqrtDetSigma + d*log(2*pi)/2;
                
            %% Univariate Gaussian model
            case 'Gaussian'

                % Parameters
                mu  = m.theta(1);
                tau = m.theta(2);

                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*tau) + (Y(I, m.Ivar) - mu).^2 / 2 / tau;                
                
            %% Poisson distribution
            case 'Poisson'
                % Parameters
                lambda = m.theta(1);
                
                % Negative log-likelihood
                subL(I) = lambda - Y(I, m.Ivar)*log(lambda) + gammaln(Y(I,m.Ivar)+1);
                
            %% geometric distribution
            case 'geometric'
                % Parameters
                theta = m.theta(1);
                
                % Negative log-likelihood
                subL(I) = -Y(I, m.Ivar)*log(1 - theta) - log(theta);
                
            %% Inverse Gaussian
            case 'invGaussian'
                % Parameters
                mu = m.theta(1);
                lambda = m.theta(2);
                
                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*lambda) + (1/2)*log(Y(I, m.Ivar).^3) - 1/lambda/mu + Y(I, m.Ivar)/2/mu^2/lambda + (1./Y(I, m.Ivar))/2/lambda;
                
            %% Gaussian linear regression
            case 'linreg'
                
                % Parameters
                tau = m.theta(1);
                b0  = m.theta(2);
                b   = m.theta(3:end);
                
                X   = Y(:,mm.ModelTypes{i}.CovIx);
                I   = I & ~any(isnan(X),2);         % missing covs
                mu  = b0 + X(I,:)*b;
                
                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*tau) + (Y(I, m.Ivar) - mu).^2 / 2 / tau;                
                
            %% Logistic regression
            case 'logreg'
                
                % Constraints
                lowerBnd = log(eps); 
                upperBnd = -lowerBnd;     
                yLims = [eps, 1-eps];                
                
                % Parameters
                b0  = m.theta(1);
                b   = m.theta(2:end);       
                
                X   = Y(:,mm.ModelTypes{i}.CovIx);
                I   = I & ~any(isnan(X),2);         % missing covs
                mu  = constrain(b0 + X(I,:)*b, lowerBnd, upperBnd);                
                
                phat=1./(1 + exp(-mu));
                if any(phat < yLims(1) | yLims(2) < phat)
                    phat = max(min(phat,yLims(2)),yLims(1));
                end
                
                % Negative log-likelihood
                subL(I) = -(Y(I, m.Ivar).*log(phat) + (1.0-Y(I, m.Ivar)).*log(1.0-phat));                
                
            otherwise
                error('Model not available');
        end
        
        L(:,k) = L(:,k) + subL;
    end
    
    % Finally, weight by mixing proportions
    L(:,k) = L(:,k) - log(a(k));
end

end