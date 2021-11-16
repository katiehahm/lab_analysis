%% CODE INCOMPLETE
%% TODO: implement better numerical stability for high-dimensional VMF
function f = besselratio(d, x)

%% Special cases [obtained from Mathematica]
if(d == 3)
    f = -(1/x) + coth(x);
elseif (d == 5)
    f = -(3/x) + x/(-1 + x*coth(x));
elseif (d == 7)
    f = -(5/x) - x/3 + x^3/(9 + 3*x^2 - 9*x*coth(x));
elseif (d == 9)
    f = (-5*x*(21 + 2*x^2)*cosh(x) + (105 + 45*x^2 + x^4)*sinh(x))/(x*(x*(15 + x^2)*cosh(x) - 3*(5 + 2*x^2)*sinh(x)));
else
    % Use lower bound in Theorem 5 [1]
    nu = d/2;
    alpha = 2;
    lam = nu + (alpha-1)/2;
    del = (nu - 1/2) + lam/2/sqrt(lam^2 + x^2);
    f = x / (del + sqrt(del^2 + x^2));
end 



end