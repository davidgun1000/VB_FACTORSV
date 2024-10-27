function [deriv]=deriv_psi_to_gamYJ(theta,gam_YJ)

    id_neg=theta<0;
    theta_bar(id_neg,1)=(1-theta(id_neg,1));
    deriv_num(id_neg,1)=((2-gam_YJ(id_neg,1)).*((theta_bar(id_neg,1)).^(2-gam_YJ(id_neg,1))).*log(theta_bar(id_neg,1)))-(((theta_bar(id_neg,1)).^(2-gam_YJ(id_neg,1)))+1);
    deriv_den(id_neg,1)=((2-gam_YJ(id_neg,1)).^2);
    deriv(id_neg,1)=deriv_num(id_neg,1)./deriv_den(id_neg,1);
    
    id_pos=theta>0;
    deriv_num(id_pos,1)=(gam_YJ(id_pos,1)).*((1+theta(id_pos,1)).^(gam_YJ(id_pos,1))).*log(1+theta(id_pos,1))-((1+theta(id_pos,1)).^(gam_YJ(id_pos,1)))+1;
    deriv_den(id_pos,1)=(gam_YJ(id_pos,1).^2);
    deriv(id_pos,1)=deriv_num(id_pos,1)./deriv_den(id_pos,1); 
    
%     if theta<0
%        theta_bar=(1-theta);
%        deriv_num=((2-gam_YJ).*((theta_bar).^(2-gam_YJ)).*log(theta_bar))-(((theta_bar).^(2-gam_YJ))+1);
%        deriv_den=((2-gam_YJ).^2);
%        deriv=deriv_num./deriv_den;
%     else
%        deriv_num=(gam_YJ).*((1+theta).^(gam_YJ)).*log(1+theta)-((1+theta).^(gam_YJ))+1;
%        deriv_den=(gam_YJ.^2);
%        deriv=deriv_num./deriv_den; 
%     end




end