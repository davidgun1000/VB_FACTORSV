function sirhat=backward_simulation_SV_noncentered(particles,w,indx,phi,sig,kapha)
%keyboard
% BACKTRACING (Kitagawa 1996, Andrieu et al, 2010)
T=size(w,1);
N=size(w,2);
outndx=NaN(1,T);
outndx(T)=randsample(N,1,true,w(T,:));
sirhat(:,T)=particles(:,outndx(T),T);
for t=T-1:-1:1,
    log_weight_backward=log(w(t,:))-0.5*log(2*pi)-0.5.*((sirhat(:,t+1)-phi.*(particles(:,:,t))).^2);
    w_backward=exp(log_weight_backward-max(log_weight_backward));
    w_backward=w_backward./sum(w_backward);
    indx(t,1)=find(rand(1) < cumsum(w_backward),1,'first');
    sirhat(:,t)=particles(:,indx(t,1),t);
end
%outndx(t)=indx(t+1,outndx(t+1));    
%log_weight_ancestral=log(w(t-1,:))-0.5*log(2*pi)-0.5*log(tau)-0.5.*(1/tau).*((particles(1,1,t)-mu-phi.*(particles(:,:,t-1)-mu)).^2);
%w_ancestral=exp(log_weight_ancestral-max(log_weight_ancestral));
%w_ancestral=w_ancestral./sum(w_ancestral);
%indx(t,1)=find(rand(1) < cumsum(w_ancestral),1,'first');