% Negative log-likelihood for the VMF distribution
function negl = vmflike(mu, kappa, x)
            
d = length(mu);
if(d == 3)
    logbesseli = log(sinh(kappa)) + 0.5*log(2) - 0.5*log(kappa) - 0.5*log(pi);
else
    % TODO: improve numerical accuracy
    logbesseli = log(besseli(d/2-1,kappa));
end

negl = sum( -( (d/2-1)*log(kappa) - d/2*log(2*pi) - logbesseli ) - kappa*x*mu );

end