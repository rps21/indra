function [logpdf] = laplacelogpdf(x, mu, b)
%LAPLACELOGPDF Calculate log-pdf (unnormalized) for laplace distribution
%
%   [logpdf] = laplacelogpdf(x, mu, b)
%   
%   where: laplace(x|mu,b) = 1/(2*b)*exp( -abs(x-mu)/b )
%    and   log(laplace(x|mu,b) ~ -abs(x-mu)/b
%
%   this function is vectorized

x = 10.^x;
mu = 10.^mu;
b = 10.^b;

logpdf = -abs(x-mu)./b;
