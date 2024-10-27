function [out]=grad_gamYJ_to_gamYJtransform(gam_YJ_transform)

    %out=(2*exp(gam_YJ_transform))./((1+exp(gam_YJ_transform)).^2);
    a=0.001;
    b=1.999;
    out1=(b-a).*exp(-gam_YJ_transform);
    out2=(1+exp(-gam_YJ_transform)).^2;
    out=out1./out2;
end