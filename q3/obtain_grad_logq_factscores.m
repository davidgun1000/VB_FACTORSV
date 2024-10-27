function [grad_logq]=obtain_grad_logq_factscores(theta,psi,mu_VB_factscores,d_VB_factscores,gam_YJ_factscores)

    grad_logq_temp1=(second_derivative_theta(theta,gam_YJ_factscores)./deriv_psi_to_theta(theta,gam_YJ_factscores));
    grad_logq_temp2=(-deriv_psi_to_theta(theta,gam_YJ_factscores)).*((psi-mu_VB_factscores)./(d_VB_factscores.^2));
    grad_logq=grad_logq_temp1+grad_logq_temp2;



end