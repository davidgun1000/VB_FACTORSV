function [particles,w,indx,sir_llh,u_res]=csmc_mu_factorSV_corr_noncentered(y,B,ft,phi,sig,kapha,N,T,ctraj,u1,num,dim)

if num<=dim
   t=1;
   particles=zeros(1,N,T);
   w=zeros(T,N);
   indx=zeros(T,N);
   particles(1,1,:)=ctraj;
   particles(:,2:N,1)=sqrt(1/(1-phi^2))*u1(2:end,t)';
   logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)+kapha))-0.5.*(1./exp(sig.*particles(1,:,t)+kapha)).*((y(1,t)-B(1,:)*ft(:,t)).^2);
   w(t,:)=exp(logw-max(logw));
   sir_llh=log(mean(w(t,:)))+max(logw);
   w(t,:)=w(t,:)./sum(w(t,:)); 
   for t=2:T
       u=sort(rand(1,N-1));            
       [indx(t,:),u_res(:,t)]=rs_multinomial_cond_corr2(particles(:,1:N,t-1),w(t-1,:),u);
       particles(:,2:N,t)=phi*(particles(:,indx(t,2:N),t-1))+u1(2:end,t)';
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
    particles(1,1,:)=ctraj;
    particles(:,2:N,1)=sqrt(1/(1-phi^2))*u1(2:end,t)';
    logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)))-0.5.*exp(-(sig.*particles(1,:,t))).*(ft(num-dim,t).^2);
    w(t,:)=exp(logw-max(logw));
    sir_llh=log(mean(w(t,:)))+max(logw);
    w(t,:)=w(t,:)./sum(w(t,:)); 
    for t=2:T
        u=sort(rand(1,N-1));            
       [indx(t,:),u_res(:,t)]=rs_multinomial_cond_corr2(particles(:,1:N,t-1),w(t-1,:),u);
        particles(:,2:N,t)=phi*(particles(:,indx(t,2:N),t-1))+u1(2:end,t)';
        logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)))-0.5.*exp(-(sig.*particles(1,:,t))).*(ft(num-dim,t).^2);
        w(t,:)=exp(logw-max(logw)); 
        sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
        w(t,:)=w(t,:)./sum(w(t,:));
    end
end
end