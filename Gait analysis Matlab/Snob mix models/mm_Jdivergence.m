% BETA - work in progress
function Jdiv = mm_Jdivergence(mm, wModels)

% If model matrix wasn't passed
if(~exist('wModels','var'))
    wModels = 1:mm.nModelTypes;
end

M = length(wModels);
K = mm.nClasses;   % number of mixtures
Jdiv = cell(M,1);

%% Compute J divergence 
% For each model...
for i = wModels
    
    Jdiv{i} = zeros(K,K);
    
    % Compute J divergence matrix between all the classes
    for j = 1:(K-1)
        for k = (j+1):K
            
            switch mm.ModelTypes{i}.type 
                % Exponential distribution
                case 'exp'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*(-2 + a/b + b/a);   % J divergence
                    
                % gamma distribution
                case 'gamma'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    mu1 = a(1); v1 = a(2);
                    mu2 = b(1); v2 = b(2);  
                    
                    f = gammaln(v1) - gammaln(v2) + v1*(psi(v1) - 1) + v2*(log(mu2*v1) - log(mu1*v2) - psi(v1) + mu1/mu2);
                    f = f + gammaln(v2) - gammaln(v1) + v2*(psi(v2) - 1) + v1*(log(mu1*v2) - log(mu2*v1) - psi(v2) + mu2/mu1);
                    f = f / 2;
                
                % Geometric distribution
                case 'geometric'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*-(((a - b)*(log(1 - a) - log(1 - b)))/(a*b));
                    
                % Poisson distribution    
                case 'Poisson'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*(a - b)*log(a/b);
                    
                % Normal distribution
                case 'Gaussian'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    mu1 = a(1); v1 = a(2);
                    mu2 = b(1); v2 = b(2);
                    
                    f = (mu1-mu2)^2/2*(1/v1+1/v2) + 1/2*(v1/v2 + v2/v1 - 2);
                
                % Inverse Gaussian distribution
                case 'invGaussian'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    mu1 = a(1); lam1 = a(2);
                    mu2 = b(1); lam2 = b(2);  
                    
                    f = 1/lam2*(lam1/2 + 1/2/mu1 + mu1/2/mu2^2 - 1/mu2) - 1/2;
                    f = f + 1/lam1*(lam2/2 + 1/2/mu2 + mu2/2/mu1^2 - 1/mu1) - 1/2;
                    f = f / 2;
                
                % Multivariate Gaussian distribution
                case 'mvg'
                    d = mm.ModelTypes{i}.nDim;  % d-variate Gaussian
                    theta1 = mm.class{j}.model{i}.theta;
                    theta2 = mm.class{k}.model{i}.theta; % Model 2 parameters                    
                    mu1 = theta1(1:d);
                    mu2 = theta2(1:d);
                    Sigma1 = reshape(theta1(d+1:end),d,d);
                    Sigma2 = reshape(theta2(d+1:end),d,d);
                    
                    f = mvgkl(mu1, Sigma1, mu2, Sigma2);
                    f = f + mvgkl(mu2, Sigma2, mu1, Sigma1);
                    f = f / 2;
                    
                % Single Factor Analysis model
                case 'sfa'
                    d = mm.ModelTypes{i}.nDim;  % d-variate Gaussian
                    theta1 = mm.class{j}.model{i}.theta; % Model 1 parameters                           
                    theta2 = mm.class{k}.model{i}.theta; % Model 2 parameters                           
                    
                    % extract parameters
                    mu1 = theta1(1:d);
                    mu2 = theta2(1:d); 
                    s1 = theta1(d+1:2*d);
                    s2 = theta2(d+1:2*d);
                    a1 = theta1((2*d+1):end);
                    a2 = theta2((2*d+1):end);   % TODO: collapsed factor
                    
                    % compute SFA covariance
                    Sigma1 = a1*a1' + diag(s1.^2);
                    Sigma2 = a2*a2' + diag(s2.^2);
                    
                    f = mvgkl(mu1, Sigma1, mu2, Sigma2);
                    f = f + mvgkl(mu2, Sigma2, mu1, Sigma1);
                    f = f / 2;                    
                                                    
                % Laplace distribution
                case 'Laplace'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters                    
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    mu1 = a(1); lam1 = a(2);
                    mu2 = b(1); lam2 = b(2);     
                    
                    absmu = abs(mu1-mu2);
                    f = (lam1*exp(-absmu/lam1) + lam2*(-log(lam1)+log(lam2) - 1) + absmu)/lam2;
                    f = f + (lam2*exp(-absmu/lam2) + lam1*(-log(lam2)+log(lam1) - 1) + absmu)/lam1;
                    f = f / 2;
                    
                % Multinomial distribution
                case 'multi'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters                    
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    
                    f = sum(a .* log(a./b));
                    f = f + sum(b .* log(b./a));
                    f = f / 2;
                    
                % Weibull distribution
                case 'weibull'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    lam1 = a(1); k1 = a(2);
                    lam2 = b(1); k2 = b(2);
                    
                    g = -psi(1);    % Euler–Mascheroni constant
                    f = log(k1) - k1*log(lam1) - log(k2) + k2*log(lam2) + (k1-k2)*(log(lam1) - g/k1) + (lam1/lam2)^k2*gamma(k2/k1+1) - 1;
                    f = f + log(k2) - k2*log(lam2) - log(k1) + k1*log(lam1) + (k2-k1)*(log(lam2) - g/k2) + (lam2/lam1)^k1*gamma(k1/k2+1) - 1;
                    f = f / 2;
                               
                otherwise
                    error('Model not available');                    
            end
            
            % Save result to matrix
            Jdiv{i}(j,k) = f;
            Jdiv{i}(k,j) = f;
        end
    end
end


end