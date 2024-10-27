function [deriv]=deriv_theta_to_gamYJ(theta,psi,gam_YJ)

      id_neg=theta<0 & psi<0;
      deriv(id_neg,1)=-deriv_theta_to_psi(psi(id_neg,1),gam_YJ(id_neg,1)).*deriv_psi_to_gamYJ(theta(id_neg,1),gam_YJ(id_neg,1));

      id_pos=theta>0 & psi>0;
      deriv(id_pos,1)=-deriv_theta_to_psi(psi(id_pos,1),gam_YJ(id_pos,1)).*deriv_psi_to_gamYJ(theta(id_pos,1),gam_YJ(id_pos,1));  
%     if theta<0 & psi<0
%        deriv=-deriv_theta_to_psi(psi,gam_YJ).*deriv_psi_to_gamYJ(theta,gam_YJ);
%     end
%     
%     if theta>0 & psi>0
%        deriv=-deriv_theta_to_psi(psi,gam_YJ).*deriv_psi_to_gamYJ(theta,gam_YJ);        
%     end

end