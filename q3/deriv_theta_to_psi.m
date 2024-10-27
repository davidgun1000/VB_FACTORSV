function [deriv]=deriv_theta_to_psi(psi,gam_YJ)

   id_neg=psi<0;
   deriv(id_neg,1)=((1-psi(id_neg,1).*(2-gam_YJ(id_neg,1))).^((gam_YJ(id_neg,1)-1)./(2-gam_YJ(id_neg,1))));
   id_pos=psi>0;
   deriv(id_pos,1)=(1+psi(id_pos,1).*gam_YJ(id_pos,1)).^((1-gam_YJ(id_pos,1))./gam_YJ(id_pos,1));  


   %if theta<0
   %   deriv=((1-psi.*(2-gam_YJ)).^((gam_YJ-1)./(2-gam_YJ)));       
   %else
   %   deriv=(1+psi.*gam_YJ).^((1-gam_YJ)./gam_YJ);  
   %end


end



