function [particles,w,indx,sir_llh]=smc_noncentered(y,phi,sig,kapha,N,T,u1,u1_res)

     t=1; 
     particles=zeros(1,N,T); 
     w=zeros(T,N);
     indx=zeros(T,N); 
     particles(:,1:N,t)=sqrt(1/(1-phi^2))*u1(:,t)';
     logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)+kapha))-0.5.*(1./exp(sig.*particles(1,:,t)+kapha)).*((y(1,t)).^2);
     w(t,:)=exp(logw-max(logw));
     sir_llh=log(mean(w(t,:)))+max(logw);
     w(t,:)=w(t,:)./sum(w(t,:));
     for t=2:T
         [indx(t,:)]=rs_multinomial_corr2_PMMH(particles(:,1:N,t-1),w(t-1,:),u1_res(:,t)');
         particles(:,1:N,t)=phi*(particles(:,indx(t,:),t-1))+u1(:,t)';
         logw=-0.5*log(2*pi)-0.5.*log(exp(sig.*particles(1,:,t)+kapha))-0.5.*(1./exp(sig.*particles(1,:,t)+kapha)).*((y(1,t)).^2);
         w(t,:)=exp(logw-max(logw));
         sir_llh=sir_llh+log(mean(w(t,:)))+max(logw);
         w(t,:)=w(t,:)./sum(w(t,:));
     end

end

