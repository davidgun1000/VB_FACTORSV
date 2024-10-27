function [result]=compute_woodbury_pred(beta_loading,idio_h,fact_lambda)

D_t_inv=diag([exp(-fact_lambda)]);
U_t_inv=diag([exp(-idio_h)]);

temp=U_t_inv*beta_loading;
result=U_t_inv-temp*((D_t_inv+beta_loading'*U_t_inv*beta_loading)\(temp'));
end