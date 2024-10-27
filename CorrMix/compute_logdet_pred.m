function [result]=compute_logdet_pred(beta_loading,idio_h,fact_lambda)

Dt_inv=diag([exp(-fact_lambda)]);
Ut_inv=diag([exp(-idio_h)]);
temp=Dt_inv+beta_loading'*Ut_inv*beta_loading;
Dt=diag([exp(fact_lambda)]);
Ut=diag([exp(idio_h)]);
result=logdet(temp,'chol')+logdet(Dt,'chol')+logdet(Ut,'chol');
end