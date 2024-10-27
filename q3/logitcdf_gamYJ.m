function [y]=logitcdf_gamYJ(x)
% input x in (-inf,inf)
%global PhiMax
%y=1./(1+exp(-x)); % transform in (0,1)
%y=PhiMax*y; % transform in (0,phimax)
   %y=(2.*exp(x))./(1+exp(x));
    a=0.001;
    b=1.999;
    y=a+(b-a).*(1./(1+exp(-x)));
end