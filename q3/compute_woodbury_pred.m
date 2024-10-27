function [result]=compute_woodbury_pred(beta_loading,idio_h,fact_lambda)

D_t_inv=diag([exp(-fact_lambda)]);
U_t_inv=diag([exp(-idio_h)]);

temp=U_t_inv*beta_loading;
[chol_temp,flag]=chol(inv((D_t_inv+beta_loading'*U_t_inv*beta_loading)));
if flag==0
    result=U_t_inv-temp*((D_t_inv+beta_loading'*U_t_inv*beta_loading)\(temp'));
else
    result=NaN;
end