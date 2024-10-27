function [log_weight]=compute_IW_factorSV(y,theta_G_idio,theta_L_idio,mu_G_idio,C_G_idio,mu_L_idio,C_L_idio,num_G_idio,num_L_idio,s_G_idio,s_L_idio,...
    theta_G_fact,theta_L_fact,mu_G_fact,C_G_fact,mu_L_fact,C_L_fact,num_G_fact,num_L_fact,s_G_fact,s_L_fact,...
    theta_beta,psi_beta,s_VB_beta,epsilon_VB_beta,mu_VB_beta,B_VB_beta,d_VB_beta,gam_YJ_beta,...
    theta_factscores,...
    dim_y,num_fact,num_param_factscores,num_param_betaloading,prior,num_factor_VB,log_weight_prev)


    T=length(y(1,:)');
    
    for i=1:dim_y
        param_kapha_idio(i,1)=theta_G_idio{i,1}(1,1);
        param_phi_idio(i,1)=exp(theta_G_idio{i,1}(2,1))/(1+exp(theta_G_idio{i,1}(2,1)));
        param_sig_idio(i,1)=log(exp(theta_G_idio{i,1}(3,1))+1);
        ctraj_idio(i,:)=theta_L_idio{i,1}(1:T,1)';
        param_psi_idio(i,1)=theta_G_idio{i,1}(2,1);
        param_alpha_idio(i,1)=theta_G_idio{i,1}(3,1);
    end

    for i=1:num_fact
        param_phi_factor(i,1)=exp(theta_G_fact{i,1}(1,1))/(1+exp(theta_G_fact{i,1}(1,1)));
        param_sig_factor(i,1)=log(exp(theta_G_fact{i,1}(2,1))+1);
        ctraj_factor(i,:)=theta_L_fact{i,1}(1:T,1)';
        param_psi_factor(i,1)=theta_G_fact{i,1}(1,1);
        param_alpha_factor(i,1)=theta_G_fact{i,1}(2,1);
    end

    beta_loading_temp=theta_beta;
    combine_beta=[];
    
    for i=1:num_fact
        fact_score(i,:)=theta_factscores{i,1};
        %psi_fact_score(i,:)=psi_theta_factscores{i,1};
        beta_loading(:,i)=[zeros(i-1,1);beta_loading_temp(((i-1)*dim_y-sum(0:(i-2)))+1:i*dim_y-sum(0:(i-1)),1)];
        beta_loading(i,i)=exp(beta_loading(i,i));
        combine_beta=[combine_beta;beta_loading(i:end,i)];        
    end
    
    %------------------
    %compute logq
    log_variational_idio=0;
    for i=1:dim_y
        log_q_theta_G_idio=-num_G_idio/2*log(2*pi)+sum(log(diag(C_G_idio{i,1})))-0.5.*(s_G_idio{i,1}'*s_G_idio{i,1});
        log_q_theta_L_idio=-num_L_idio/2*log(2*pi)+sum(log(diag(C_L_idio{i,1})))-0.5.*(s_L_idio{i,1}'*s_L_idio{i,1});
        log_q_theta_idio=log_q_theta_G_idio+log_q_theta_L_idio;
        log_variational_idio=log_variational_idio+log_q_theta_idio;
    end

    log_variational_fact=0;
    for i=1:num_fact
        log_q_theta_G_fact=-num_G_fact/2*log(2*pi)+sum(log(diag(C_G_fact{i,1})))-0.5.*(s_G_fact{i,1}'*s_G_fact{i,1});
        log_q_theta_L_fact=-num_L_fact/2*log(2*pi)+sum(log(diag(C_L_fact{i,1})))-0.5.*(s_L_fact{i,1}'*s_L_fact{i,1});
        log_q_theta_fact=log_q_theta_G_fact+log_q_theta_L_fact;
        log_variational_fact=log_variational_fact+log_q_theta_fact;
    end

%     log_variational_factscores=0;
%     for i=1:num_fact
%        log_q_factscores_term1=-(T/2)*log(2*pi)-0.5.*sum(log(d_VB_factscores{i,1}.^2))-0.5.*sum(((psi_theta_factscores{i,1}-mu_VB_factscores{i,1}).^2)./(d_VB_factscores{i,1}.^2));
%        log_q_factscores_term2=sum(log(deriv_psi_to_theta(fact_score(i,:)',gam_YJ_VB_factscores{i,1})));
%        log_q_factscores=log_q_factscores_term1+log_q_factscores_term2;
%        log_variational_factscores=log_variational_factscores+log_q_factscores;
%         
%     end
    
    
    temp1=sum(log(d_VB_beta.^2))+logdet(eye(num_factor_VB)+B_VB_beta'*(diag([1./(d_VB_beta.^2)]))*B_VB_beta);
    temp_wood_beta=compute_woodbury(B_VB_beta,d_VB_beta);
    temp_compute_beta=(temp_wood_beta*(B_VB_beta*s_VB_beta+d_VB_beta.*epsilon_VB_beta));
    temp2_betaloading=sum(log(deriv_psi_to_theta(theta_beta,gam_YJ_beta)));
    log_variational_betaloading=-(num_param_betaloading/2)*log(2*pi)-0.5*temp1-0.5*(B_VB_beta*s_VB_beta+d_VB_beta.*epsilon_VB_beta)'*temp_compute_beta+temp2_betaloading;

    log_variational=log_variational_idio+log_variational_fact+log_variational_betaloading;
    
    %----------------------
    %compute log post
    
    for t=1:T
        logdet_covmat=compute_logdet_pred(beta_loading,param_sig_idio.*ctraj_idio(:,t)+param_kapha_idio,param_sig_factor.*ctraj_factor(:,t));
        inv_covmat=compute_woodbury_pred(beta_loading,param_sig_idio.*ctraj_idio(:,t)+param_kapha_idio,param_sig_factor.*ctraj_factor(:,t));
        temp_log_dens(t,1)=-(dim_y/2)*log(2*pi)-0.5*logdet_covmat-0.5*y(:,t)'*inv_covmat*y(:,t);
    end
    log_posterior1=sum(temp_log_dens);
    %if sum(isnan(log_posterior1))>0
    %   log_weight=log_weight_prev; 
    %else
    
    log_posterior4=0;
    for i=1:dim_y
        log_posterior4=log_posterior4-0.5.*log(2*pi)-0.5.*log(1./(1-param_phi_idio(i,1)^2))-0.5.*(1-param_phi_idio(i,1).^2).*(ctraj_idio(i,1).^2)-...
            (T-1)/2*log(2*pi)-0.5.*sum((ctraj_idio(i,2:T)-param_phi_idio(i,1).*ctraj_idio(i,1:T-1)).^2);        
    end
    
    log_posterior3=0;
    for i=1:num_fact
        log_posterior3=log_posterior3-0.5.*log(2*pi)-0.5.*log(1./(1-param_phi_factor(i,1).^2))-0.5.*(1-param_phi_factor(i,1).^2).*(ctraj_factor(i,1).^2)-...
            (T-1)/2*log(2*pi)-0.5.*sum((ctraj_factor(i,2:T)-param_phi_factor(i,1).*ctraj_factor(i,1:T-1)).^2);
    end
    
    log_posterior=log_posterior3+log_posterior4+log_posterior1;
    
    diag_betaloadingprior=0;
    for i=1:num_fact
        diag_betaloadingprior=diag_betaloadingprior+sum(log(beta_loading(i,i)));
    end
    
    log_prior_param_psi_idio=sum((prior.a0-1).*log((1+param_phi_idio)./2)+(prior.b0-1).*log((1-param_phi_idio)./2))+sum(param_psi_idio-log((1+exp(param_psi_idio)).^2));
    log_prior_param_psi_factor=sum((prior.a0-1).*log((1+param_phi_factor)./2)+(prior.b0-1).*log((1-param_phi_factor)./2))+sum(param_psi_factor-log((1+exp(param_psi_factor)).^2));
    log_prior_param_alpha_idio=sum(log(2)-log(pi)-log(1+param_sig_idio.^2))+sum(param_alpha_idio-log(1+exp(param_alpha_idio)));
    log_prior_param_alpha_factor=sum(log(2)-log(pi)-log(1+param_sig_factor.^2))+sum(param_alpha_factor-log(1+exp(param_alpha_factor)));
    log_prior_kapha_idio=sum(log(normpdf(param_kapha_idio,0,sqrt(prior.hp_sig2))));
    log_prior_betaloading=sum(log(normpdf(combine_beta,0,1)))+diag_betaloadingprior;
    
    
    
    log_prior=log_prior_param_psi_idio+log_prior_param_psi_factor+log_prior_param_alpha_idio+log_prior_param_alpha_factor+log_prior_kapha_idio+log_prior_betaloading;
    log_posterior_prior=log_posterior+log_prior;
    
    log_weight=log_posterior_prior-log_variational;
    %end
end