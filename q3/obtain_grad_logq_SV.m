function [grad_logq]=obtain_grad_logq_SV(d_L,D_L,F_L,C_L,C_G,theta_L,s_G,s_L,D2_star,indx_CL,func_v_eye_numL)
    
    temp2=(theta_L(indx_CL(:,1),1)-d_L(indx_CL(:,1),1)).*s_L(indx_CL(:,2),1);    
    temp=(func_v_eye_numL-D2_star.*temp2);    
    grad_logq_thetaG_term1=F_L'*temp;
    grad_logq_thetaG_term2=-C_G*s_G;
    grad_logq_thetaG_term3=-D_L'*s_L;
    grad_logq_thetaG=grad_logq_thetaG_term1+grad_logq_thetaG_term2+grad_logq_thetaG_term3;
    grad_logq_thetaL=-C_L*s_L;
    grad_logq=[grad_logq_thetaG;grad_logq_thetaL];    

end