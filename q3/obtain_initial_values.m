%obtain initial values

function [mu_VB_idio,mu_VB_fact,mu_VB_factscores,mu_VB_beta]=obtain_initial_values(phi_idio,mu_idio,tau_idio,phi_fact,tau_fact,dim_y,num_fact,length_time)

t=1;
h=sqrt(tau_idio/(1-phi_idio^2))*randn(dim_y,1)+mu_idio;
ctraj_idio_init(:,t)=h;
hf=sqrt(tau_fact/(1-phi_fact^2))*randn(num_fact,1);
ctraj_factor_init(:,t)=hf;

for t=2:length_time
    h=mu_idio+phi_idio*(h-mu_idio)+sqrt(tau_idio)*randn(dim_y,1);
    ctraj_idio_init(:,t)=h;    
    hf=phi_fact*hf+sqrt(tau_fact)*randn(num_fact,1);
    ctraj_factor_init(:,t)=hf;    
end

for i=1:dim_y
    mu_VB_idio{i,1}=[ctraj_idio_init(i,:)';mu_idio;log(phi_idio/(1-phi_idio));log(tau_idio)];    
end

for i=1:num_fact
    mu_VB_fact{i,1}=[ctraj_factor_init(i,:)';log(phi_fact/(1-phi_fact));log(tau_fact)];
end

f_init=exp(ctraj_factor_init./2).*randn(num_fact,length_time);

for i=1:num_fact
    mu_VB_factscores{i,1}=f_init(i,:)';
end

mu_VB_beta=[];
for i=1:num_fact
    mu_VB_beta=[mu_VB_beta;rand(dim_y-(i-1),1)];
end

%mu_VB_betafact=[];
%for i=1:num_fact
%    mu_VB_betafact=[mu_VB_betafact;f_init(i,:)'];
%end
%for i=1:num_fact
%    mu_VB_betafact=[mu_VB_betafact;rand(dim_y-(i-1),1)]; 
%end

%mu_VB_betafact=[mu_VB_betafact;0.5;rand(dim_y-1,1);0;rand(dim_y-2,1);0.2;rand(dim_y-3,1)]; 



%mu_VB_betafact=[f_init(1,:)';f_init(2,:)';f_init(3,:)';f_init(4,:)';log(beta_loading_true(1,1));
%    beta_loading_true(2:end,1);log(beta_loading_true(2,2));beta_loading_true(3:end,2);
%    log(beta_loading_true(3,3));beta_loading_true(4:end,3);
%    log(beta_loading_true(4,4));beta_loading_true(5:end,4)];
        


end