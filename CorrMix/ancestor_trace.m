function sirhat=ancestor_trace(particles,w,indx)
%particles: 1 by N by T
%w: T by N
%indx: T by N
% BACKTRACING (Kitagawa 1996, Andrieu et al, 2010)
T=size(w,1);
N=size(w,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1
    outndx(t)=indx(t+1,outndx(t+1));
    sirhat(:,t)=particles(:,outndx(t),t);
end