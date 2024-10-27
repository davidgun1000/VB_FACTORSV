function [logpdf]=logcauchy(tau)

logpdf=-0.5*log(tau)-log(pi*(1+tau));