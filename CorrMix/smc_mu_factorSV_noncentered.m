function [particles,w,indx,sir_llh]=smc_mu_factorSV_noncentered(y,B,ft,phi,sig,kapha,N,T,num,dim)

if num<=dim
   t=1;
   particles=zeros(1,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(:,1:N,1)=sqrt(1/(1-phi^2))*randn(1,N);
   logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)+kapha))-0.5.*(1./exp(sig.*particles(1,:,t)+kapha)).*((y(1,t)-B(1,:)*ft(:,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:)); 
   for t=2:T
       indx(t,:)=rs_multinomial(w(t-1,:));
       particles(:,1:N,t)=phi*(particles(:,indx(t,1:N),t-1))+randn(1,N);
       logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)+kapha))-0.5.*(1./exp(sig.*particles(1,:,t)+kapha)).*((y(1,t)-B(1,:)*ft(:,t)).^2);
       w(t,:)=exp(logw-max(logw)); 
       sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
       w(t,:)=w(t,:)./sum(w(t,:));
   end
 
else
    t=1;
    particles=zeros(1,N,T);
    w=zeros(T,N);
    indx=zeros(T,N);
    particles(:,1:N,1)=sqrt(1/(1-phi^2))*randn(1,N);
    logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)))-0.5.*exp(-(sig.*particles(1,:,t))).*(ft(num-dim,t).^2);
    w(t,:)=exp(logw-max(logw));
    sir_llh=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 
    for t=2:T
        indx(t,:)=rs_multinomial(w(t-1,:));
        particles(:,1:N,t)=phi*(particles(:,indx(t,1:N),t-1))+randn(1,N);
        logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)))-0.5.*exp(-(sig.*particles(1,:,t))).*(ft(num-dim,t).^2);
        w(t,:)=exp(logw-max(logw)); 
        sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
        w(t,:)=w(t,:)./sum(w(t,:));
    end
end
end