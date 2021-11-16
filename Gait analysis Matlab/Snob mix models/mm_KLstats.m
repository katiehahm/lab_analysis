%MM_KLSTATS    Print matrices of KL divergences for each combination of classes.
%  MM_KLSTATS(.) prints a (d+1)x(d+1) matrix of Kullback-Leilber (KL) divergences
%  for each attribute in the data set, where d is the number of classes. 
%  Entry (i+1,j+1) in the KL matrix corresponds to the KL divergence from class i to class j
%  for an attribute. 
%  The class labelled "Pop" represents a mixture model with a single class
%  only - that is, a mixture model where all data belongs to one class.
%
%  The function does nothing if the passed mixture model mm only contains a
%  single class.
%  
%  The input arguments are:
%   mm    - structure respresenting the complete mixture model
%   data  - data set that was used to fit the mixture model
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function mm_KLstats(mm, data)

nClasses = mm.nClasses;
if(nClasses > 1)
    % Fit a mixture model with one class for the entire sample
    Pop = snob(data, mm.opts.ModelList, 'k', 1, 'fixedstructure', true, 'display', false);

    % For each attribute
    nModels = mm.nModelTypes;
    for i = 1:nModels
        klmat = zeros(nClasses+1, nClasses+1);
        
        % Compute KL matrix between class j and k for attr i 
        for j = 0:nClasses
            for k = 0:nClasses
                if(j ~= k)
                    klmat(j+1,k+1) = ComputeKL(data, Pop, mm, i, j, k);
                end
            end
        end

        %% Print KL matrix for attribute i
        Ivar  = mm.ModelTypes{i}.Ivar;
        D     = length(Ivar);
        fprintf('** ');
        if(D > 1) 
            fprintf('['); 
        end
        for j = 1:D
            fprintf('%s', mm.opts.VarNames{Ivar(j)})
            if(j<D)
                fprintf(' ');
            end
        end
        if(D > 1) 
            fprintf(']'); 
        end
        fprintf(' ~ %s\n', mm.ModelTypes{i}.Description);
        printMat(klmat);
        fprintf('\n');
    end

end


end

%% Compute KL between two models
function kl = ComputeKL(data, Pop, mm, attr, classSrc, classTgt)

% Attribute type
type = mm.ModelTypes{attr}.type;

% Parameters of source and target model
thetaSrc = Pop.class{1}.model{attr}.theta;
if(classSrc > 0)
    thetaSrc = mm.class{classSrc}.model{attr}.theta;
end
thetaTgt = Pop.class{1}.model{attr}.theta;
if(classTgt > 0)
    thetaTgt = mm.class{classTgt}.model{attr}.theta;
end

% Compute KL
switch type
    case 'beta'
        a0 = thetaSrc(1); b0 = thetaSrc(2);
        a1 = thetaTgt(1); b1 = thetaTgt(2);

        kl = betaln(a1,b1)-betaln(a0,b0) + (a0-a1)*psi(a0) + (b0-b1)*psi(b0) + (a1-a0+b1-b0)*psi(a0+b0);

    case 'cfixexp'
        c = mm.ModelTypes{attr}.c;
        alpha = thetaSrc; beta = thetaTgt;
        kl = (exp(-c./alpha)-1).*(1-alpha./beta+log(alpha./beta));

    case 'cfixweibull'
        c = mm.ModelTypes{attr}.c; 
        lambda1 = thetaSrc(1); k1 = thetaSrc(2);
        lambda2 = thetaTgt(1); k2 = thetaTgt(2);        
        kl = klwbltypeI(c, k1, lambda1, k2, lambda2);

    case 'exp'
        p0 = thetaSrc; p1 = thetaTgt;
        kl = -1 + p0/p1 - log(p0) + log(p1);
    
    case 'gamma'
        mu0 = thetaSrc(1); phi0 = thetaSrc(2);
        k0 = phi0; theta0 = mu0 / phi0;
        mu1 = thetaTgt(1); phi1 = thetaTgt(2);
        k1 = phi1; theta1 = mu1 / phi1;

        kl = -k0 + (k0*theta0)/theta1 + k1*log(theta1/theta0) - gammaln(k0) + gammaln(k1) + (k0 - k1)*psi(k0);
    
    case 'Gaussian'
        m0 = thetaSrc(1); v0 = thetaSrc(2);
        m1 = thetaTgt(1); v1 = thetaTgt(2);
        kl = 0.5*(v0/v1 + (m1 - m0)^2/v1 + log(v1/v0) - 1);

    case 'geometric'
        p = thetaSrc; q = thetaTgt;
        kl = log(p/q) - (1 - p)/p*log((1 - q)/(1 - p));

    case 'invGaussian'
        mu1 = thetaSrc(1); lambda1 = thetaSrc(2);
        mu2 = thetaTgt(1); lambda2 = thetaTgt(2); 

        kl = 0.5*log(lambda2/lambda1) + 1/lambda2*(lambda1/2 + 1/2/mu1 + mu1/2/mu2^2 - 1/mu2) - 0.5;

    case 'Laplace'
        m0 = thetaSrc(1); v0 = thetaSrc(2);
        m1 = thetaTgt(1); v1 = thetaTgt(2);        
        kl = 1/v1*(abs(m0-m1) + v1*log(2*v1) + v0*exp(-abs(m0-m1)/v0)) - (1 + log(2*v0));
           
    case 'linreg'
        X  = data(:,mm.ModelTypes{attr}.CovIx);

        v0 = thetaSrc(1);
        b0 = thetaSrc(2);
        b  = thetaSrc(3:end);
        m0 = b0 + X*b;

        v1 = thetaTgt(1);
        b0 = thetaTgt(2);
        b  = thetaTgt(3:end);        
        m1 = b0 + X*b;

        kl = mean( 0.5*(v0/v1 + (m1 - m0).^2/v1 + log(v1/v0) - 1) );

    case 'logreg'
        X  = data(:,mm.ModelTypes{attr}.CovIx);

        b0 = thetaSrc(1);
        b  = thetaSrc(2:end);
        m0 = logsig(b0 + X*b);
        p0 = 1./(1 + exp(-m0));
        p0 = max(eps, p0); p0 = min(1-eps, p0);

        b0 = thetaTgt(1);
        b  = thetaTgt(2:end);        
        m1 = b0 + X*b;
        p1 = 1./(1 + exp(-m1));
        p1 = max(eps, p1); p1 = min(1-eps, p1);

        kl = mean( p0.*log(p0./p1) + (1-p0).*log((1-p0)./(1-p1)) ); 

    case 'multi'
        kl = sum(thetaSrc .* log(thetaSrc ./ thetaTgt));

    case 'negb'
        Nmax = 1e6;      
        mu = thetaSrc(1); phi = thetaSrc(2);
        r0 = phi; p0 = 1-mu/(mu+phi);
        mu = thetaTgt(1); phi = thetaTgt(2);        
        r1 = phi; p1 = 1-mu/(mu+phi);
        kl = -(nbinlike([r0,p0], nbinrnd(r0,p0,Nmax,1)) - nbinlike([r1,p1], nbinrnd(r0,p0,Nmax,1))) ./ Nmax;        
    
    case 'mvg'
        d = mm.ModelTypes{attr}.nDim;

        m0 = thetaSrc(1:d);
        S0 = reshape(thetaSrc(d+1:end),d,d);
        m1 = thetaTgt(1:d);
        S1 = reshape(thetaTgt(d+1:end),d,d);        

        kl = mvgkl(m0, S0, m1, S1) / d; % per-dimension

    case 'Poisson'
        p0 = thetaSrc; p1 = thetaTgt;
        kl = -p0 + p1 + p0*log(p0) - p0*log(p1);

    case 'sfa'
        d = mm.ModelTypes{attr}.nDim;

        mu0 = thetaSrc(1:d);      
        sigma = thetaSrc(d+1:2*d); 
        a_sfa = thetaSrc(2*d+1:end); 
        Sigma0 = a_sfa*a_sfa' + diag(sigma.^2);

        mu1 = thetaTgt(1:d);      
        sigma = thetaTgt(d+1:2*d); 
        a_sfa = thetaTgt(2*d+1:end); 
        Sigma1 = a_sfa*a_sfa' + diag(sigma.^2);     

        kl = mvgkl(mu0, Sigma0, mu1, Sigma1) / d; % per-dimension

    case 'vmf'
        kappa0 = thetaSrc(1);
        mu0 = thetaSrc(2:end);
        kappa1 = thetaTgt(1);
        mu1 = thetaTgt(2:end);    

        Nmax = 1e5;             
        kl = -(vmflike(mu0, kappa0, vmfrnd(kappa0,mu0,Nmax)) - vmflike(mu1,kappa1, vmfrnd(kappa0,mu0,Nmax))) ./ Nmax;

    case 'weibull'
        lambda1 = thetaSrc(1); k1 = thetaSrc(2);
        lambda2 = thetaTgt(1); k2 = thetaTgt(2);        
        
        g = -psi(1);
        kl = log(k1/lambda1^k1) - log(k2/lambda2^k2) + (k1 - k2)*(log(lambda1) -  g/k1) + (lambda1/lambda2)^k2*gamma(k2/k1 + 1) - 1;

    otherwise
        error('KL divergence not implemented this model type yet');
end


end

%% Print the KL matrix
function printMat(x)

n = size(x,1);

fprintf('%15s','');
for i = 1:n
    str = 'Pop.';
    if(i>1)
        str = ['Class ',num2str(i-1)];
    end
    fprintf('%15s', str);
end
fprintf('\n');

for i = 1:n
    str = 'Pop.';
    if(i>1)
        str = ['Class ',num2str(i-1)];
    end
    fprintf('%15s', str);
    for j = 1:n
        if(i~=j)
            if(x(i,j) > 1000)
                fprintf('%15s', '>1000');
            else
                fprintf('%15.4f', x(i,j));
            end
        else
            fprintf('%15s', '-');
        end
    end
    fprintf('\n');
end


end

% KL Divergence of Weibull Type I with Censoring
function klest = klwbltypeI(c, k, lambda, k1, lambda1)

EulerGamma = -psi(1);
p = 1 - exp(-(c/lambda)^k);
p = max(p, eps);
p = min(p, 1-eps);
Ei = ExpIntegralEi(log(1-p));

A = k*log(c/lambda1) + 1/(1-p)*(EulerGamma - Ei + k*log(lambda1/lambda));
B = 1 - (1+EulerGamma)/(1-p) - log(k) + k*log(lambda/c) + log(k1) + 1/(1-p)*(Ei + log(k/k1) + (lambda/lambda1)^k1*(gamma(1+k1/k) - MathGamma(-log(1-p),1+k1/k))) + (c/lambda1)^k1;
klest = (1-p)/k*(k1*A + k*B);

end

function f = ExpIntegralEi(x)

f = real(-expint(-x));

end

function f = MathGamma(z, a)

f = gammainc(z, a, 'upper') * gamma(a);

end
