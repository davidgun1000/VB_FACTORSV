function [u1_used]=obtain_random_numbers(phi,sig,kapha,ctraj,N,T)

for t=1:T
    if t==1
        u1(1,t)=sqrt((1-phi^2))*(ctraj(1,1));
    else
        u1(1,t)=(ctraj(1,t)-phi*(ctraj(1,t-1)));
    end
end

u1_used=[u1;randn(N-1,T)]; 
end
