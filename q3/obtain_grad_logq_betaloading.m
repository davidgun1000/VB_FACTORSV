function [grad_logq]=obtain_grad_logq_betaloading(theta,psi,mu_VB_beta,B_VB_beta,d_VB_beta,gam_YJ,temp_wood_beta)

    grad_logq_temp1=(second_derivative_theta(theta,gam_YJ)./deriv_psi_to_theta(theta,gam_YJ));
    grad_logq_temp2=(-deriv_psi_to_theta(theta,gam_YJ)).*(temp_wood_beta*(psi-mu_VB_beta));
    grad_logq=grad_logq_temp1+grad_logq_temp2;



end