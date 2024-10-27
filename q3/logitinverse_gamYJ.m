function y=logitinverse_gamYJ(x)
% input x in (0,phimax)

%x=x/PhiMax; % transform in (0,1)
%y=log((x./2)./(1-(x./2))); % transform to (-inf,inf)
a=0.001;
b=1.999;
y=log((x-a)./(b-a))-log(1-((x-a)./(b-a)));

end
