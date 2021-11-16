function [lse,sm] = logsumexp(x)
%LOGSUMEXP  Log-sum-exp function.
%    lse = LOGSUMEXP(x) returns the log-sum-exp function evaluated at 
%    the vector x, defined by lse = log(sum(exp(x)).
%    [lse,sm] = LOGSUMEXP(x) also returns the softmax function evaluated
%    at x, defined by sm = exp(x)/sum(exp(x)).
%
%    NOTE: this is a vectorized version of:
%    https://au.mathworks.com/matlabcentral/fileexchange/84892-logsumexp-softmax
% 
%    Reference:
%    P. Blanchard, D. J. Higham, and N. J. Higham.  
%    Accurately computing the log-sum-exp and softmax functions. 
%    IMA J. Numer. Anal., Advance access, 2020.


xmax = max(x,[],2);
e = exp(x - xmax);
s = sum(e, 2) - 1;           % not ideal ...
lse = xmax + log1p(s);

if nargout > 1
    sm = e/(1+s);
end 

end