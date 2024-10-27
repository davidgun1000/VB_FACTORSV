function [deriv]=deriv_psi_to_theta(theta,gam_YJ)

   id_neg=theta<0;
   deriv(id_neg,1)=((1-theta(id_neg,1)).^(1-gam_YJ(id_neg,1))); 
   id_pos=theta>0;
   deriv(id_pos,1)=(1+theta(id_pos,1)).^(gam_YJ(id_pos,1)-1);   

%    if theta<0
%       deriv=((1-theta)^(1-gam_YJ)); 
%    else
%       deriv=(1+theta)^(gam_YJ-1);        
%    end

end