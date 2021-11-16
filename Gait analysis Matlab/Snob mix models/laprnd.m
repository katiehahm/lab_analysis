% Generate random numbers from the Laplace distribution
function x = laprnd(mu, b, n)

U = rand(n,1) - 0.5;
x = mu - b*sign(U).*log(1-2*abs(U));

end