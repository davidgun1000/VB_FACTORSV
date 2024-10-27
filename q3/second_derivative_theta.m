function [deriv]=second_derivative_theta(theta,gam_YJ)

    id_neg=theta<0;
    deriv(id_neg,1)=(gam_YJ(id_neg,1)-1).*((1-theta(id_neg,1)).^(-gam_YJ(id_neg,1)));
    id_pos=theta>0;
    deriv(id_pos,1)=(gam_YJ(id_pos,1)-1).*((1+theta(id_pos,1)).^(gam_YJ(id_pos,1)-2));


end