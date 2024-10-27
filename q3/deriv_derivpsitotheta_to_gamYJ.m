function [deriv]=deriv_derivpsitotheta_to_gamYJ(theta,gam_YJ)

    id_neg=theta<0;
    theta_bar(id_neg,1)=1-theta(id_neg,1);
    deriv(id_neg,1)=-((theta_bar(id_neg,1)).^(1-gam_YJ(id_neg,1))).*log(theta_bar(id_neg,1));    

    id_pos=theta>0;
    deriv(id_pos,1)=((theta(id_pos,1)+1).^(gam_YJ(id_pos,1)-1)).*log(theta(id_pos,1)+1); 
    %if theta<0
    %    theta_bar=1-theta;
    %    deriv=-((theta_bar).^(1-gam_YJ)).*log(theta_bar);       
    %else
    %    deriv=((theta+1).^(gam_YJ-1)).*log(theta+1);                
    %end
    
end