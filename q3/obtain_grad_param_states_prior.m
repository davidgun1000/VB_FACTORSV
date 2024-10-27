function [grad_param_states]=obtain_grad_param_states_prior(theta,y,num_param,num_states,prior)
    
   param_kapha=theta(1,1);
   param_phi=exp(theta(2,1))/(1+exp(theta(2,1)));
   param_sig=log(exp(theta(3,1))+1);
   
   T=length(y);
   grad_kapha=0.5.*(sum((y.^2).*exp(-param_sig.*theta(num_param+1:num_param+T,1)-param_kapha))-T)-param_kapha/prior.hp_sig2;
   %grad_psi=(sum((theta(num_param+2:num_param+T,1)-param_phi.*(theta(num_param+1:num_param+T-1,1))).*(theta(num_param+1:num_param+T-1,1)))+(theta(num_param+1,1)^2)*param_phi-...
   %    (param_phi/(1-param_phi^2)))*(param_phi*(1-param_phi))-theta(2,1)/prior.hp_sig2;
   
   grad_psi_prior_term1=(((prior.a0-1)./(1+param_phi))-((prior.b0-1)./(1-param_phi))).*(exp(theta(2,1))./((1+exp(theta(2,1))).^2));
   grad_psi_prior_term2=((1-exp(theta(2,1)))./(1+exp(theta(2,1))));
   grad_psi_prior=grad_psi_prior_term1+grad_psi_prior_term2;
   
   grad_psi=(sum((theta(num_param+2:num_param+T,1)-param_phi.*(theta(num_param+1:num_param+T-1,1))).*(theta(num_param+1:num_param+T-1,1)))+(theta(num_param+1,1)^2)*param_phi-...
       (param_phi/(1-param_phi^2)))*(param_phi*(1-param_phi))+grad_psi_prior;

   grad_alpha_prior_term1=((-2.*param_sig)./(1+param_sig.^2)).*(exp(theta(3,1))./(1+exp(theta(3,1))));
   grad_alpha_prior_term2=1-(exp(theta(3,1))./(1+exp(theta(3,1))));
   grad_alpha_prior=grad_alpha_prior_term1+grad_alpha_prior_term2;
   
   grad_alpha=0.5.*(sum(theta(num_param+1:num_param+T,1).*(y.^2).*exp(-param_sig.*theta(num_param+1:num_param+T,1)-param_kapha)-theta(num_param+1:num_param+T,1)).*(1-exp(-param_sig)))-...
       grad_alpha_prior;
   
   %grad_alpha=0.5.*(sum(theta(num_param+1:num_param+T,1).*(y.^2).*exp(-param_sig.*theta(num_param+1:num_param+T,1)-param_kapha)-theta(num_param+1:num_param+T,1)).*(1-exp(-param_sig)))-...
   %    theta(3,1)/prior.hp_sig2;
   
   grad_states=zeros(T,1);
   grad_states(1,1)=(param_sig/2)*((y(1,1)^2)*exp(-param_sig*theta(num_param+1,1)-param_kapha)-1)+param_phi*(theta(num_param+2,1)-param_phi*theta(num_param+1,1))-...
       theta(num_param+1,1)*((1-param_phi)^2);
   
   grad_states(T,1)=(param_sig/2)*((y(T,1)^2)*exp(-param_sig*theta(num_param+T,1)-param_kapha)-1)-(theta(num_param+T,1)-param_phi*theta(num_param+T-1,1));
   grad_states(2:T-1,1)=(param_sig/2).*((y(2:T-1,1).^2).*exp(-param_sig*theta(num_param+2:num_param+T-1,1)-param_kapha)-1)+param_phi.*(theta(num_param+3:num_param+T,1)-param_phi.*theta(num_param+2:num_param+T-1,1))-...
       (theta(num_param+2:num_param+T-1,1)-param_phi.*theta(num_param+1:num_param+T-2,1));
  
   grad_param=[grad_kapha;grad_psi;grad_alpha];
   grad_param_states=[grad_param;grad_states];
end

% grad_psi=(-(param_phi/((1-param_phi^2)))+(param_phi/param_tau)*(beta(1,1)-param_mu)^2+...
%          (1/param_tau).*sum((beta(2:T,1)-param_mu-param_phi.*(beta(1:T-1,1)-param_mu)).*(beta(1:T-1,1)-param_mu)))*...
%          (exp(beta(num_states+2,1))/((1+exp(beta(num_states+2,1)))^2))+...
%          1-((2*exp(beta(num_states+2,1)))/(1+exp(beta(num_states+2,1))))+(exp(beta(num_states+2,1))/((1+exp(beta(num_states+2,1)))^2))*...
%          (((prior.a0-1)/(1+param_phi))-((prior.b0-1)/(1-param_phi)));



