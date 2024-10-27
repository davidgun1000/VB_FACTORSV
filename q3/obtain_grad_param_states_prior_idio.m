function [grad_param_states]=obtain_grad_param_states_prior_idio(theta_G,theta_L,y,theta_factscores,theta_beta,num_param_idio,num_states_idio,numer,prior,dim_y,num_fact)
    
    param_kapha=theta_G(1,1);
    param_phi=exp(theta_G(2,1))/(1+exp(theta_G(2,1)));
    param_sig=log(exp(theta_G(3,1))+1);
    param_psi=theta_G(2,1);
    param_alpha=theta_G(3,1);
    for i=1:num_fact
        fact(i,:)=theta_factscores{i,1};
        beta_loading(:,i)=[zeros(i-1,1);theta_beta(((i-1)*dim_y-sum(0:(i-2)))+1:i*dim_y-sum(0:(i-1)),1)];
        beta_loading(i,i)=exp(beta_loading(i,i));
        
    end
    
    
    T=length(y(1,:)');
    
    grad_kapha=0.5.*(sum(((y(numer,:)-beta_loading(numer,:)*fact)'.^2).*exp(-param_sig.*theta_L(1:T,1)-param_kapha))-T)-param_kapha./prior.hp_sig2;
    
    grad_psi_prior_term1=(((prior.a0-1)./(1+param_phi))-((prior.b0-1)./(1-param_phi))).*(exp(theta_G(2,1))./((1+exp(theta_G(2,1))).^2));
    grad_psi_prior_term2=((1-exp(theta_G(2,1)))./(1+exp(theta_G(2,1))));
    grad_psi_prior=grad_psi_prior_term1+grad_psi_prior_term2;
    
    grad_psi=(sum((theta_L(2:T,1)-param_phi.*(theta_L(1:T-1,1))).*(theta_L(1:T-1,1)))+(theta_L(1,1)^2).*param_phi-...
        (param_phi./(1-param_phi^2)))*(param_phi*(1-param_phi))+grad_psi_prior;
    
    grad_alpha_prior_term1=((-2.*param_sig)./(1+param_sig.^2)).*(exp(theta_G(3,1))./(1+exp(theta_G(3,1))));
    grad_alpha_prior_term2=1-(exp(theta_G(3,1))./(1+exp(theta_G(3,1))));
    grad_alpha_prior=grad_alpha_prior_term1+grad_alpha_prior_term2;
    
    grad_alpha=0.5.*(sum(theta_L(1:T,1).*((y(numer,:)-beta_loading(numer,:)*fact)'.^2).*exp(-param_sig.*theta_L(1:T,1)-param_kapha)-theta_L(1:T,1)).*(1-exp(-param_sig)))+...
        grad_alpha_prior;

    grad_states=zeros(T,1);
    grad_states(1,1)=(param_sig./2).*(((y(numer,1)-beta_loading(numer,:)*fact(:,1)).^2).*exp(-param_sig*theta_L(1,1)-param_kapha)-1)+param_phi.*(theta_L(2,1)-param_phi.*theta_L(1,1))-...
        theta_L(1,1)*((1-param_phi)^2);    
    grad_states(T,1)=(param_sig./2).*(((y(numer,T)-beta_loading(numer,:)*fact(:,T)).^2).*exp(-param_sig*theta_L(T,1)-param_kapha)-1)-(theta_L(T,1)-param_phi.*theta_L(T-1,1));
    grad_states(2:T-1,1)=(param_sig/2).*(((y(numer,2:T-1)-beta_loading(numer,:)*fact(:,2:T-1))'.^2).*exp(-param_sig.*theta_L(2:T-1,1)-param_kapha)-1)+param_phi.*(theta_L(3:T,1)-param_phi.*theta_L(2:T-1,1))-...
        (theta_L(2:T-1,1)-param_phi.*theta_L(1:T-2,1));
   
    grad_param=[grad_kapha;grad_psi;grad_alpha];
    grad_param_states=[grad_param;grad_states];
end








