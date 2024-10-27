load('DATA_SP100.mat'); %load the data
T_est=1000; % length of time series
y=y_used_demean(:,1:T_est); % the dataset
num_L_idio=length(y(1,:)');
num_L_fact=length(y(1,:)');
num_G_idio=3;
num_G_fact=2;
dim_y=size(y,1); % the dimension of y.
num_fact=1; % the number of latent factors
T=length(y(1,:)');

indx_CG_idio=[1,1;2,1;3,1;2,2;3,2;3,3];
indx_CG_fact=[1,1;2,1;2,2];
[indx_CL_idio]=obtain_index_C2(num_L_idio);
[indx_CL_fact]=obtain_index_C2(num_L_fact);
[indx_f_F_idio]=obtain_index_f_F(num_L_idio);
[indx_f_F_fact]=obtain_index_f_F(num_L_fact);

indx_diag_CG_idio(:,1)=(1:num_G_idio)';
indx_diag_CG_idio(:,2)=(1:num_G_idio)';

indx_diag_CG_fact(:,1)=(1:num_G_fact)';
indx_diag_CG_fact(:,2)=(1:num_G_fact)';

indx_diag_CL_idio(:,1)=(1:num_L_idio)';
indx_diag_CL_idio(:,2)=(1:num_L_idio)';

indx_diag_CL_fact(:,1)=(1:num_L_fact)';
indx_diag_CL_fact(:,2)=(1:num_L_fact)';

func_v_eye_numL_idio=func_v(eye(num_L_idio));
func_v_eye_numL_idio(indx_f_F_idio,1)=0.001;
id=func_v_eye_numL_idio~=0.001;
func_v_eye_numL_idio=func_v_eye_numL_idio(id,1);

func_v_eye_numL_fact=func_v(eye(num_L_fact));
func_v_eye_numL_fact(indx_f_F_fact,1)=0.001;
id=func_v_eye_numL_fact~=0.001;
func_v_eye_numL_fact=func_v_eye_numL_fact(id,1);

A_temp=[0;1];
repmat_01_idio=repmat(A_temp,2*num_L_idio,1);
repmat_01_idio=repmat_01_idio(1:2*num_L_idio-1,1);
repmat_01_fact=repmat(A_temp,2*num_L_fact,1);
repmat_01_fact=repmat_01_fact(1:2*num_L_fact-1,1);

num_factor_VB=4; % the number of VB factors. 
num_param_factscores=T;
num_param_betaloading=sum(1:dim_y)-sum(1:dim_y-num_fact); % the number of parameters in the factor loading. 


%priors
prior.a0=20;
prior.b0=1.5;
prior.hp_sig2=10;

const_A_idio=2*num_L_idio-1;
const_A_fact=2*num_L_fact-1;

%ADAM hyperparameters
adapt_tau_1=0.9;
adapt_tau_2=0.99;
adapt_epsilon=10^-8;
adapt_alpha=0.001;

%initial values

for i=1:dim_y
    mu_G_idio{i,1}=[0;4;-3];
    C_G_idio{i,1}=4*speye(num_G_idio);
    C_G_star_idio{i,1}=4*speye(num_G_idio);
    C_L_star_idio{i,1}=speye(num_L_idio);
    C_L_idio{i,1}=speye(num_L_idio);
    d_L_idio{i,1}=zeros(num_L_idio,1);
    D_L_idio{i,1}=zeros(num_L_idio,num_G_idio);
    f_L_idio{i,1}=0.1*rand(const_A_idio,1);
    F_L_idio{i,1}=0.1*rand(const_A_idio,num_G_idio);
    m_t_idio{i,1}=0;
    v_t_idio{i,1}=0;  
end 
length_mu_G_idio=length(mu_G_idio{1,1});
length_C_G_idio=length(func_v(C_G_idio{1,1}));
length_d_L_idio=length(d_L_idio{1,1});
length_D_L_idio=length(vec(D_L_idio{1,1}));
length_f_L_idio=length(f_L_idio{1,1});
length_F_L_idio=length(vec(F_L_idio{1,1}));


for i=1:num_fact
    mu_G_fact{i,1}=[4;-2];
    C_G_fact{i,1}=4*speye(num_G_fact);
    C_G_star_fact{i,1}=4*speye(num_G_fact);
    C_L_star_fact{i,1}=speye(num_L_fact);
    C_L_fact{i,1}=speye(num_L_fact);
    if i==1
        d_L_fact{i,1}=0*ones(num_L_fact,1);
    else
        d_L_fact{i,1}=-10*ones(num_L_fact,1);
    end
    D_L_fact{i,1}=zeros(num_L_fact,num_G_fact);
    f_L_fact{i,1}=0.1*rand(const_A_fact,1);
    F_L_fact{i,1}=0.1*rand(const_A_fact,num_G_fact);
    m_t_fact{i,1}=0;
    v_t_fact{i,1}=0;
    
end

length_mu_G_fact=length(mu_G_fact{1,1});
length_C_G_fact=length(func_v(C_G_fact{1,1}));
length_d_L_fact=length(d_L_fact{1,1});
length_D_L_fact=length(vec(D_L_fact{1,1}));
length_f_L_fact=length(f_L_fact{1,1});
length_F_L_fact=length(vec(F_L_fact{1,1}));

mu_VB_beta=[];
for i=1:num_fact
    mu_VB_beta=[mu_VB_beta;zeros(dim_y-(i-1),1)];
end
d_VB_beta=0.01*ones(num_param_betaloading,1);
for i=1:num_factor_VB
    B_VB_beta(:,i)=[zeros(i-1,1);0.01*rand(num_param_betaloading-i+1,1)];
end
gam_YJ_beta=ones(num_param_betaloading,1);
gam_YJ_beta_transform=logitinverse_gamYJ(gam_YJ_beta);

length_mu_VB_beta=length(mu_VB_beta);
length_B_VB_beta=length(func_v(B_VB_beta));
length_d_VB_beta=length(d_VB_beta);
length_gam_YJ_beta_transform=length(gam_YJ_beta_transform);

m_t_beta=0;
v_t_beta=0;

D_square_inv_factscores=sparse(num_param_factscores,num_param_factscores);
indx_diag_d_factscores(:,1)=(1:num_param_factscores)';
indx_diag_d_factscores(:,2)=(1:num_param_factscores)'; 

D_square_inv_beta=sparse(num_param_betaloading,num_param_betaloading);
indx_diag_d_beta(:,1)=(1:num_param_betaloading)';

n_iter=50000; %the number of VB iterations

s_idio=cell(dim_y,1);
theta_idio=cell(dim_y,1);
 
s_fact=cell(dim_y,1);
theta_fact=cell(dim_y,1);

log_weight_prev=0;
for i=1:n_iter
    i
    %tic
    
    %generating Monte Carlo samples
    
    s_VB_beta=randn(num_factor_VB,1);
    epsilon_VB_beta=randn(num_param_betaloading,1);
    psi_beta=mu_VB_beta+B_VB_beta*s_VB_beta+d_VB_beta.*epsilon_VB_beta;
    gam_YJ_beta=logitcdf_gamYJ(gam_YJ_beta_transform);
    theta_beta=YJ_psi_to_theta(psi_beta,gam_YJ_beta);
    
    mu_VB_beta_prev=mu_VB_beta;
    B_VB_beta_prev=B_VB_beta;
    d_VB_beta_prev=d_VB_beta;
    gam_YJ_beta_prev=gam_YJ_beta;

    for j=1:dim_y
        s_idio{j,1}=randn(num_L_idio+num_G_idio,1);
        s_G_idio{j,1}=s_idio{j,1}(1:num_G_idio,1);
        s_L_idio{j,1}=s_idio{j,1}(num_G_idio+1:end,1);
        var1_idio{j,1}=((C_G_idio{j,1}')\s_G_idio{j,1});
        theta_G_idio{j,1}=mu_G_idio{j,1}+var1_idio{j,1};
        C_G_idio_prev{j,1}=C_G_idio{j,1};
        v_C_L_star_idio{j,1}=f_L_idio{j,1}+F_L_idio{j,1}*theta_G_idio{j,1};
        id=v_C_L_star_idio{j,1}~=0;
        v_C_L_temp_idio{j,1}=v_C_L_star_idio{j,1}(id,1);
         
        C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_CL_idio(:,1)'],[indx_CL_idio(:,2)']))=v_C_L_temp_idio{j,1};
        C_L_idio{j,1}=C_L_star_idio{j,1};
        temp_C_L_star_idio{j,1}=exp(C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)'])));
        C_L_idio{j,1}(sub2ind(size(C_L_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)']))=temp_C_L_star_idio{j,1}';
        C_L_idio_prev{j,1}=C_L_idio{j,1};
        var2_idio{j,1}=((C_L_idio{j,1}')\s_L_idio{j,1});
        var3_idio{j,1}=((C_L_idio{j,1}')\D_L_idio{j,1});
        var3_var1_idio{j,1}=var3_idio{j,1}*var1_idio{j,1};
        mu_L_idio{j,1}=d_L_idio{j,1}-var3_var1_idio{j,1};        
        theta_L_idio_temp=mu_L_idio{j,1}+var2_idio{j,1};
        id_idio=abs(theta_L_idio_temp)>=80;
        if sum(id_idio)>0
           theta_L_idio{j,1}=theta_L_idio{j,1};
        else
           theta_L_idio{j,1}=theta_L_idio_temp;
        end
    end
    
    for j=1:num_fact
        s_fact{j,1}=randn(num_L_fact+num_G_fact,1);
        s_G_fact{j,1}=s_fact{j,1}(1:num_G_fact,1);
        s_L_fact{j,1}=s_fact{j,1}(num_G_fact+1:end,1);
        var1_fact{j,1}=((C_G_fact{j,1}')\s_G_fact{j,1});
        theta_G_fact{j,1}=mu_G_fact{j,1}+var1_fact{j,1};
        C_G_fact_prev{j,1}=C_G_fact{j,1};
        v_C_L_star_fact{j,1}=f_L_fact{j,1}+F_L_fact{j,1}*theta_G_fact{j,1};
        id=v_C_L_star_fact{j,1}~=0;
        v_C_L_temp_fact{j,1}=v_C_L_star_fact{j,1}(id,1);
         
        C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_CL_fact(:,1)'],[indx_CL_fact(:,2)']))=v_C_L_temp_fact{j,1};
        C_L_fact{j,1}=C_L_star_fact{j,1};
        temp_C_L_star_fact{j,1}=exp(C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)'])));
        C_L_fact{j,1}(sub2ind(size(C_L_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)']))=temp_C_L_star_fact{j,1}';
        C_L_fact_prev{j,1}=C_L_fact{j,1};
        var2_fact{j,1}=((C_L_fact{j,1}')\s_L_fact{j,1});
        var3_fact{j,1}=((C_L_fact{j,1}')\D_L_fact{j,1});
        var3_var1_fact{j,1}=var3_fact{j,1}*var1_fact{j,1};
        mu_L_fact{j,1}=d_L_fact{j,1}-var3_var1_fact{j,1};         
        theta_L_fact_temp=mu_L_fact{j,1}+var2_fact{j,1};

        id_fact=abs(theta_L_fact_temp)>=80;
        if sum(id_fact)>0
           theta_L_fact{j,1}=theta_L_fact{j,1}; 
        else
           theta_L_fact{j,1}=theta_L_fact_temp;            
        end
    end
    
    ctraj_idio=[];
    ctraj_fact=[];
    for j=1:dim_y
        ctraj_idio=[ctraj_idio;theta_L_idio{j,1}']; 
    end
     
    for j=1:num_fact
        ctraj_fact=[ctraj_fact;theta_L_fact{j,1}'];
    end
    
    beta_loading_temp = theta_beta;
    
    for ss=1:num_fact
        beta_loading(:,ss)=[zeros(ss-1,1);beta_loading_temp(((ss-1)*dim_y-sum(0:(ss-2)))+1:ss*dim_y-sum(0:(ss-1)),1)];
        beta_loading(ss,ss)=exp(beta_loading(ss,ss));
    end
    
    for ss=1:dim_y
        temp_kapha_idio(ss,1)=theta_G_idio{ss,1}(1,1);
        temp_sig_idio(ss,1)=log(exp(theta_G_idio{ss,1}(3,1))+1); 
        
    end
    
    for ss=1:num_fact
        temp_sig_fact(ss,1) = log(exp(theta_G_fact{ss,1}(2,1))+1); 
    end
        for t=1:T
        ctraj_idio_temp(:,t)=ctraj_idio(1:dim_y,t);
        ctraj_fact_temp(:,t)=ctraj_fact(1:num_fact,t);
        
        Vt_inv=diag(1./(exp(temp_sig_idio(1:dim_y,1).*ctraj_idio_temp(:,t)+temp_kapha_idio(1:dim_y,1))));
        Dt_inv=diag(exp(-(temp_sig_fact(:,1).*ctraj_fact_temp(:,t))));
        [chol_var_ft,flag]=chol(inv((beta_loading'*Vt_inv*beta_loading)+Dt_inv),'lower');
        if flag==0
           mean_ft=(chol_var_ft*chol_var_ft')*beta_loading'*(Vt_inv*y(:,t));
           fact_score_temp(:,t)=mvnrnd(mean_ft,chol_var_ft*chol_var_ft');      
           id_fact=abs(fact_score_temp(:,t))>20;
        if sum(sum(id_fact))>0
           fact_score(:,t)=randn(num_fact,1); 
        else
           fact_score(:,t)=fact_score_temp(:,t);
        
        end
        
        
        end
     end
    
     for ss=1:num_fact
         theta_factscores{ss,1} = fact_score_temp(ss,:)';
     end
     %compute the lower bound
     [log_weight(i,1)]=compute_IW_factorSV(y,theta_G_idio,theta_L_idio,mu_G_idio,C_G_idio_prev,mu_L_idio,C_L_idio_prev,num_G_idio,num_L_idio,s_G_idio,s_L_idio,...
                theta_G_fact,theta_L_fact,mu_G_fact,C_G_fact_prev,mu_L_fact,C_L_fact_prev,num_G_fact,num_L_fact,s_G_fact,s_L_fact,...
                theta_beta,psi_beta,s_VB_beta,epsilon_VB_beta,mu_VB_beta_prev,B_VB_beta_prev,d_VB_beta_prev,gam_YJ_beta_prev,...
                theta_factscores,... 
                dim_y,num_fact,num_param_factscores,num_param_betaloading,prior,num_factor_VB,log_weight_prev);
            
     while (isnan(log_weight(i,1))) | isinf(log_weight(i,1)) | log_weight(i,1)<log_weight(1,1) | log_weight(i,1)>0
            s_VB_beta=randn(num_factor_VB,1);
            epsilon_VB_beta=randn(num_param_betaloading,1);
            psi_beta=mu_VB_beta+B_VB_beta*s_VB_beta+d_VB_beta.*epsilon_VB_beta;
            gam_YJ_beta=logitcdf_gamYJ(gam_YJ_beta_transform);
            theta_beta=YJ_psi_to_theta(psi_beta,gam_YJ_beta);

            mu_VB_beta_prev=mu_VB_beta;
            B_VB_beta_prev=B_VB_beta;
            d_VB_beta_prev=d_VB_beta;
            gam_YJ_beta_prev=gam_YJ_beta;

            for j=1:dim_y
                s_idio{j,1}=randn(num_L_idio+num_G_idio,1);
                s_G_idio{j,1}=s_idio{j,1}(1:num_G_idio,1);
                s_L_idio{j,1}=s_idio{j,1}(num_G_idio+1:end,1);
                var1_idio{j,1}=((C_G_idio{j,1}')\s_G_idio{j,1});
                theta_G_idio{j,1}=mu_G_idio{j,1}+var1_idio{j,1};
                C_G_idio_prev{j,1}=C_G_idio{j,1};
                v_C_L_star_idio{j,1}=f_L_idio{j,1}+F_L_idio{j,1}*theta_G_idio{j,1};
                id=v_C_L_star_idio{j,1}~=0;
                v_C_L_temp_idio{j,1}=v_C_L_star_idio{j,1}(id,1);

                C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_CL_idio(:,1)'],[indx_CL_idio(:,2)']))=v_C_L_temp_idio{j,1};
                C_L_idio{j,1}=C_L_star_idio{j,1};
                temp_C_L_star_idio{j,1}=exp(C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)'])));
                C_L_idio{j,1}(sub2ind(size(C_L_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)']))=temp_C_L_star_idio{j,1}';
                C_L_idio_prev{j,1}=C_L_idio{j,1};
                var2_idio{j,1}=((C_L_idio{j,1}')\s_L_idio{j,1});
                var3_idio{j,1}=((C_L_idio{j,1}')\D_L_idio{j,1});
                var3_var1_idio{j,1}=var3_idio{j,1}*var1_idio{j,1};
                mu_L_idio{j,1}=d_L_idio{j,1}-var3_var1_idio{j,1};        
                theta_L_idio_temp=mu_L_idio{j,1}+var2_idio{j,1};
                id_idio=abs(theta_L_idio_temp)>=80;
                if sum(id_idio)>0
                   theta_L_idio{j,1}=theta_L_idio{j,1};
                else
                   theta_L_idio{j,1}=theta_L_idio_temp;
                end
            end

            for j=1:num_fact
                s_fact{j,1}=randn(num_L_fact+num_G_fact,1);
                s_G_fact{j,1}=s_fact{j,1}(1:num_G_fact,1);
                s_L_fact{j,1}=s_fact{j,1}(num_G_fact+1:end,1);
                var1_fact{j,1}=((C_G_fact{j,1}')\s_G_fact{j,1});
                theta_G_fact{j,1}=mu_G_fact{j,1}+var1_fact{j,1};
                C_G_fact_prev{j,1}=C_G_fact{j,1};
                v_C_L_star_fact{j,1}=f_L_fact{j,1}+F_L_fact{j,1}*theta_G_fact{j,1};
                id=v_C_L_star_fact{j,1}~=0;
                v_C_L_temp_fact{j,1}=v_C_L_star_fact{j,1}(id,1);

                C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_CL_fact(:,1)'],[indx_CL_fact(:,2)']))=v_C_L_temp_fact{j,1};
                C_L_fact{j,1}=C_L_star_fact{j,1};
                temp_C_L_star_fact{j,1}=exp(C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)'])));
                C_L_fact{j,1}(sub2ind(size(C_L_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)']))=temp_C_L_star_fact{j,1}';
                C_L_fact_prev{j,1}=C_L_fact{j,1};
                var2_fact{j,1}=((C_L_fact{j,1}')\s_L_fact{j,1});
                var3_fact{j,1}=((C_L_fact{j,1}')\D_L_fact{j,1});
                var3_var1_fact{j,1}=var3_fact{j,1}*var1_fact{j,1};
                mu_L_fact{j,1}=d_L_fact{j,1}-var3_var1_fact{j,1};         
                theta_L_fact_temp=mu_L_fact{j,1}+var2_fact{j,1};

                id_fact=abs(theta_L_fact_temp)>=80;
                if sum(id_fact)>0
                   theta_L_fact{j,1}=theta_L_fact{j,1}; 
                else
                   theta_L_fact{j,1}=theta_L_fact_temp;            
                end
            end

            ctraj_idio=[];
            ctraj_fact=[];
            for j=1:dim_y
                ctraj_idio=[ctraj_idio;theta_L_idio{j,1}']; 
            end

            for j=1:num_fact
                ctraj_fact=[ctraj_fact;theta_L_fact{j,1}'];
            end

            beta_loading_temp = theta_beta;

            for ss=1:num_fact
                beta_loading(:,ss)=[zeros(ss-1,1);beta_loading_temp(((ss-1)*dim_y-sum(0:(ss-2)))+1:ss*dim_y-sum(0:(ss-1)),1)];
                beta_loading(ss,ss)=exp(beta_loading(ss,ss));
            end

            for ss=1:dim_y
                temp_kapha_idio(ss,1)=theta_G_idio{ss,1}(1,1);
                temp_sig_idio(ss,1)=log(exp(theta_G_idio{ss,1}(3,1))+1); 

            end

            for ss=1:num_fact
                temp_sig_fact(ss,1) = log(exp(theta_G_fact{ss,1}(2,1))+1); 
            end
            for t=1:T
                ctraj_idio_temp(:,t)=ctraj_idio(1:dim_y,t);
                ctraj_fact_temp(:,t)=ctraj_fact(1:num_fact,t);

                Vt_inv=diag(1./(exp(temp_sig_idio(1:dim_y,1).*ctraj_idio_temp(:,t)+temp_kapha_idio(1:dim_y,1))));
                Dt_inv=diag(exp(-(temp_sig_fact(:,1).*ctraj_fact_temp(:,t))));
                %var_ft=inv((beta_loading'*Vt_inv*beta_loading)+Dt_inv);
                [chol_var_ft,flag]=chol(inv((beta_loading'*Vt_inv*beta_loading)+Dt_inv),'lower');
                %[var_ft]=jitChol(var_ft);
                %chol_var_ft=chol(var_ft,'lower');
                if flag==0
                mean_ft=(chol_var_ft*chol_var_ft')*beta_loading'*(Vt_inv*y(:,t));
                fact_score_temp(:,t)=mvnrnd(mean_ft,chol_var_ft*chol_var_ft');      
                id_fact=abs(fact_score_temp(:,t))>20;
                if sum(sum(id_fact))>0
                   fact_score(:,t)=randn(num_fact,1); 
                else
                   fact_score(:,t)=fact_score_temp(:,t);

                end


                end
             end

             for ss=1:num_fact
                 theta_factscores{ss,1} = fact_score_temp(ss,:)';
             end
          
          
          
          
          
          
          
          
          
          
          
          [log_weight(i,1)]=compute_IW_factorSV(y,theta_G_idio,theta_L_idio,mu_G_idio,C_G_idio_prev,mu_L_idio,C_L_idio_prev,num_G_idio,num_L_idio,s_G_idio,s_L_idio,...
                theta_G_fact,theta_L_fact,mu_G_fact,C_G_fact_prev,mu_L_fact,C_L_fact_prev,num_G_fact,num_L_fact,s_G_fact,s_L_fact,...
                theta_beta,psi_beta,s_VB_beta,epsilon_VB_beta,mu_VB_beta_prev,B_VB_beta_prev,d_VB_beta_prev,gam_YJ_beta_prev,...
                theta_factscores,... 
                dim_y,num_fact,num_param_factscores,num_param_betaloading,prior,num_factor_VB,log_weight_prev);
                    
      end
     
     
     
     
     
     
     
     
     
     
     
     
    %updating the variational parameters for idiosynchratic log volatilities 
    for j=1:dim_y
         
         D1_star_idio=func_v(func_dg(C_G_idio{j,1})+ones(num_G_idio,num_G_idio)-speye(num_G_idio));         
         diag_CL_idio=diag(C_L_idio{j,1});
         diag_CL_kron_idio=kron(diag_CL_idio,[1;0]);
         diag_CL_kron_idio=diag_CL_kron_idio(1:2*num_L_idio-1,1);
         D2_star_idio=diag_CL_kron_idio+repmat_01_idio;
         [grad_param_states_idio]=obtain_grad_param_states_prior_idio(theta_G_idio{j,1},theta_L_idio{j,1},y,theta_factscores,theta_beta,num_G_idio,num_L_idio,j,prior,dim_y,num_fact);
         [grad_logq_idio]=obtain_grad_logq_SV(d_L_idio{j,1},D_L_idio{j,1},F_L_idio{j,1},C_L_idio{j,1},C_G_idio{j,1},theta_L_idio{j,1},s_G_idio{j,1},s_L_idio{j,1},...
             D2_star_idio,indx_CL_idio,func_v_eye_numL_idio);
          G_idio=grad_param_states_idio-grad_logq_idio;
          G_G_idio=G_idio(1:num_G_idio,1);
          G_L_idio=G_idio(num_G_idio+1:end,1);
         
          var4_idio=((C_L_idio{j,1})\G_L_idio);
          temp_GL_CL_idio=((G_L_idio')/(C_L_idio{j,1}'));
          temp_temp_idio=(var2_idio{j,1}(indx_CL_idio(:,1),1)-var3_var1_idio{j,1}(indx_CL_idio(:,1),1)).*(temp_GL_CL_idio(1,indx_CL_idio(:,2)))';
         
          grad_thetaL_G2_idio=-D2_star_idio.*temp_temp_idio;
          grad_mu_thetaL_G2_idio=F_L_idio{j,1}'*grad_thetaL_G2_idio;
         
          grad_LB_temp1_idio=G_G_idio+grad_mu_thetaL_G2_idio;
          grad_LB_temp2_idio=-D1_star_idio.*func_v(var1_idio{j,1}*(((G_G_idio+grad_mu_thetaL_G2_idio-(D_L_idio{j,1}')*var4_idio)')/(C_G_idio{j,1}')));
          grad_LB_temp3_idio=G_L_idio;
          grad_LB_temp4_idio=-vec(var4_idio*var1_idio{j,1}');
          grad_LB_temp5_idio=grad_thetaL_G2_idio;
          grad_LB_temp6_idio=vec(grad_thetaL_G2_idio*(theta_G_idio{j,1}'));
         
          grad_LB_idio=[grad_LB_temp1_idio;grad_LB_temp2_idio;grad_LB_temp3_idio;grad_LB_temp4_idio;grad_LB_temp5_idio;grad_LB_temp6_idio];        
          if sum(isinf(grad_LB_idio))>0 | sum(isnan(grad_LB_idio))>0
             
          else
          m_t_idio{j,1}=adapt_tau_1*m_t_idio{j,1}+(1-adapt_tau_1)*grad_LB_idio;
          v_t_idio{j,1}=adapt_tau_2*v_t_idio{j,1}+(1-adapt_tau_2)*(grad_LB_idio.^2);
          mt_hat_idio{j,1}=m_t_idio{j,1}./(1-(adapt_tau_1.^i));
          vt_hat_idio{j,1}=v_t_idio{j,1}./(1-(adapt_tau_2.^i));
         
          lambda_idio=[mu_G_idio{j,1};func_v(C_G_star_idio{j,1});d_L_idio{j,1};vec(D_L_idio{j,1});f_L_idio{j,1};vec(F_L_idio{j,1})];        
          lambda_idio=lambda_idio+adapt_alpha*(mt_hat_idio{j,1}./(sqrt(vt_hat_idio{j,1})+adapt_epsilon));
         
          mu_G_idio{j,1}=lambda_idio(1:length_mu_G_idio,1);
          temp_C_G_idio=lambda_idio(length_mu_G_idio+1:length_mu_G_idio+length_C_G_idio,1);
          d_L_idio{j,1}=lambda_idio(length_mu_G_idio+length_C_G_idio+1:length_mu_G_idio+length_C_G_idio+length_d_L_idio,1);
          D_L_idio{j,1}=inverse_vec(lambda_idio(length_mu_G_idio+length_C_G_idio+length_d_L_idio+1:length_mu_G_idio+length_C_G_idio+length_d_L_idio+length_D_L_idio,1),num_L_idio,num_G_idio);
          f_L_idio{j,1}=lambda_idio(length_mu_G_idio+length_C_G_idio+length_d_L_idio+length_D_L_idio+1:length_mu_G_idio+length_C_G_idio+length_d_L_idio+length_D_L_idio+length_f_L_idio,1);
          F_L_idio{j,1}=inverse_vec(lambda_idio(length_mu_G_idio+length_C_G_idio+length_d_L_idio+length_D_L_idio+length_f_L_idio+1:length_mu_G_idio+length_C_G_idio+length_d_L_idio+length_D_L_idio+length_f_L_idio+length_F_L_idio,1),const_A_idio,num_G_idio);
 
          C_G_star_idio{j,1}(sub2ind(size(C_G_star_idio{j,1}),[indx_CG_idio(:,1)'],[indx_CG_idio(:,2)']))=temp_C_G_idio;
          C_G_idio{j,1}=C_G_star_idio{j,1};
          temp_C_G_star_idio=exp(C_G_star_idio{j,1}(sub2ind(size(C_G_star_idio{j,1}),[indx_diag_CG_idio(:,1)'],[indx_diag_CG_idio(:,2)'])));
          C_G_idio{j,1}(sub2ind(size(C_G_idio{j,1}),[indx_diag_CG_idio(:,1)'],[indx_diag_CG_idio(:,2)']))=temp_C_G_star_idio';
         end
    end
     
    %updating the variational parameters for factor log volatilities 

    
     for j=1:num_fact
         D1_star_fact{j,1}=func_v(func_dg(C_G_fact{j,1})+ones(num_G_fact,num_G_fact)-speye(num_G_fact));
         
         diag_CL_fact{j,1}=diag(C_L_fact{j,1});
         diag_CL_kron_fact{j,1}=kron(diag_CL_fact{j,1},[1;0]);
         diag_CL_kron_fact{j,1}=diag_CL_kron_fact{j,1}(1:2*num_L_fact-1,1);
         
         D2_star_fact{j,1}=diag_CL_kron_fact{j,1}+repmat_01_fact;
         [grad_param_states_fact{j,1}]=obtain_grad_param_states_prior_fact(theta_G_fact{j,1},theta_L_fact{j,1},y,theta_factscores,theta_beta,num_G_fact,num_L_fact,j,prior,dim_y,num_fact);
         [grad_logq_fact{j,1}]=obtain_grad_logq_SV(d_L_fact{j,1},D_L_fact{j,1},F_L_fact{j,1},C_L_fact{j,1},C_G_fact{j,1},theta_L_fact{j,1},s_G_fact{j,1},s_L_fact{j,1},...
             D2_star_fact{j,1},indx_CL_fact,func_v_eye_numL_fact);
         G_fact{j,1}=grad_param_states_fact{j,1}-grad_logq_fact{j,1};
         G_G_fact{j,1}=G_fact{j,1}(1:num_G_fact,1);
         G_L_fact{j,1}=G_fact{j,1}(num_G_fact+1:end,1);
         
         var4_fact{j,1}=((C_L_fact{j,1})\G_L_fact{j,1});
         temp_GL_CL_fact{j,1}=((G_L_fact{j,1}')/(C_L_fact{j,1}'));
         temp_temp_fact{j,1}=(var2_fact{j,1}(indx_CL_fact(:,1),1)-var3_var1_fact{j,1}(indx_CL_fact(:,1),1)).*(temp_GL_CL_fact{j,1}(1,indx_CL_fact(:,2)))';
         
         grad_thetaL_G2_fact{j,1}=-D2_star_fact{j,1}.*temp_temp_fact{j,1};
         grad_mu_thetaL_G2_fact{j,1}=F_L_fact{j,1}'*grad_thetaL_G2_fact{j,1};
         
         grad_LB_temp1_fact{j,1}=G_G_fact{j,1}+grad_mu_thetaL_G2_fact{j,1};
         grad_LB_temp2_fact{j,1}=-D1_star_fact{j,1}.*func_v(var1_fact{j,1}*(((G_G_fact{j,1}+grad_mu_thetaL_G2_fact{j,1}-(D_L_fact{j,1}')*var4_fact{j,1})')/(C_G_fact{j,1}')));
         grad_LB_temp3_fact{j,1}=G_L_fact{j,1};
         grad_LB_temp4_fact{j,1}=-vec(var4_fact{j,1}*var1_fact{j,1}');
         grad_LB_temp5_fact{j,1}=grad_thetaL_G2_fact{j,1};
         grad_LB_temp6_fact{j,1}=vec(grad_thetaL_G2_fact{j,1}*(theta_G_fact{j,1}'));
         
         grad_LB_fact{j,1}=[grad_LB_temp1_fact{j,1};grad_LB_temp2_fact{j,1};grad_LB_temp3_fact{j,1};grad_LB_temp4_fact{j,1};grad_LB_temp5_fact{j,1};grad_LB_temp6_fact{j,1}];
         if sum(isinf(grad_LB_fact{j,1}))>0 | sum(isnan(grad_LB_fact{j,1}))>0
         else
         m_t_fact{j,1}=adapt_tau_1*m_t_fact{j,1}+(1-adapt_tau_1)*grad_LB_fact{j,1};
         v_t_fact{j,1}=adapt_tau_2*v_t_fact{j,1}+(1-adapt_tau_2)*(grad_LB_fact{j,1}.^2);
         mt_hat_fact{j,1}=m_t_fact{j,1}./(1-(adapt_tau_1.^i));
         vt_hat_fact{j,1}=v_t_fact{j,1}./(1-(adapt_tau_2.^i));
         
         lambda_fact{j,1}=[mu_G_fact{j,1};func_v(C_G_star_fact{j,1});d_L_fact{j,1};vec(D_L_fact{j,1});f_L_fact{j,1};vec(F_L_fact{j,1})];        
         lambda_fact{j,1}=lambda_fact{j,1}+adapt_alpha*(mt_hat_fact{j,1}./(sqrt(vt_hat_fact{j,1})+adapt_epsilon));
         
         mu_G_fact{j,1}=lambda_fact{j,1}(1:length_mu_G_fact,1);
         temp_C_G_fact{j,1}=lambda_fact{j,1}(length_mu_G_fact+1:length_mu_G_fact+length_C_G_fact,1);
         d_L_fact{j,1}=lambda_fact{j,1}(length_mu_G_fact+length_C_G_fact+1:length_mu_G_fact+length_C_G_fact+length_d_L_fact,1);
         D_L_fact{j,1}=inverse_vec(lambda_fact{j,1}(length_mu_G_fact+length_C_G_fact+length_d_L_fact+1:length_mu_G_fact+length_C_G_fact+length_d_L_fact+length_D_L_fact,1),num_L_fact,num_G_fact);
         f_L_fact{j,1}=lambda_fact{j,1}(length_mu_G_fact+length_C_G_fact+length_d_L_fact+length_D_L_fact+1:length_mu_G_fact+length_C_G_fact+length_d_L_fact+length_D_L_fact+length_f_L_fact,1);
         F_L_fact{j,1}=inverse_vec(lambda_fact{j,1}(length_mu_G_fact+length_C_G_fact+length_d_L_fact+length_D_L_fact+length_f_L_fact+1:length_mu_G_fact+length_C_G_fact+length_d_L_fact+length_D_L_fact+length_f_L_fact+length_F_L_fact,1),const_A_fact,num_G_fact);
 
         C_G_star_fact{j,1}(sub2ind(size(C_G_star_fact{j,1}),[indx_CG_fact(:,1)'],[indx_CG_fact(:,2)']))=temp_C_G_fact{j,1};
         C_G_fact{j,1}=C_G_star_fact{j,1};
         temp_C_G_star_fact{j,1}=exp(C_G_star_fact{j,1}(sub2ind(size(C_G_star_fact{j,1}),[indx_diag_CG_fact(:,1)'],[indx_diag_CG_fact(:,2)'])));
         C_G_fact{j,1}(sub2ind(size(C_G_fact{j,1}),[indx_diag_CG_fact(:,1)'],[indx_diag_CG_fact(:,2)']))=temp_C_G_star_fact{j,1}';
         end
     end
     

    
    %updating beta
    
    [grad_betaloading]=obtain_grad_factor_beta_loading_SV(theta_G_idio,theta_G_fact,theta_beta,theta_factscores,y,...
        ctraj_idio,ctraj_fact,num_param_factscores,num_param_betaloading,dim_y,num_fact);
    if sum(sum(isinf(grad_betaloading)))>0 | sum(sum(isnan(grad_betaloading)))>0
    else
    temp_wood_beta=compute_woodbury(B_VB_beta,d_VB_beta);
    [grad_logq_betaloading]=obtain_grad_logq_betaloading(theta_beta,psi_beta,mu_VB_beta,B_VB_beta,d_VB_beta,gam_YJ_beta,temp_wood_beta);
    diff_gradlogpos_gradlogq_betaloading=grad_betaloading-grad_logq_betaloading;
    
    temp_deriv_theta_to_psi_betaloading=deriv_theta_to_psi(psi_beta,gam_YJ_beta);
    grad_mu_VB_betaloading=temp_deriv_theta_to_psi_betaloading.*diff_gradlogpos_gradlogq_betaloading;
    diff_gradlogpos_gradlogq_betaloading_repmat=repmat(diff_gradlogpos_gradlogq_betaloading,num_factor_VB,1);
    s_VB_beta_kron=kron(s_VB_beta,ones(num_param_betaloading,1));
    temp_deriv_theta_to_psi_betaloading_repmat=repmat(temp_deriv_theta_to_psi_betaloading,num_factor_VB,1);
    grad_B_VB_betaloading=(temp_deriv_theta_to_psi_betaloading_repmat.*s_VB_beta_kron).*diff_gradlogpos_gradlogq_betaloading_repmat;
    grad_B_VB_betaloading=reshape(grad_B_VB_betaloading,num_param_betaloading,num_factor_VB);
    for j=2:num_factor_VB
        grad_B_VB_betaloading(1:j-1,j)=0;
    end
    grad_d_VB_betaloading=(temp_deriv_theta_to_psi_betaloading.*epsilon_VB_beta).*diff_gradlogpos_gradlogq_betaloading; 
    grad_gamYJ_betaloading_transform=(deriv_theta_to_gamYJ(theta_beta,psi_beta,gam_YJ_beta)).*diff_gradlogpos_gradlogq_betaloading.*grad_gamYJ_to_gamYJtransform(gam_YJ_beta_transform);
    
    grad_LB_beta=[grad_mu_VB_betaloading;func_v(grad_B_VB_betaloading);grad_d_VB_betaloading;grad_gamYJ_betaloading_transform];
    if sum(sum(isinf(grad_LB_beta)))>0 | sum(sum(isnan(grad_LB_beta)))>0
    else
    m_t_beta=adapt_tau_1*m_t_beta+(1-adapt_tau_1)*grad_LB_beta;
    v_t_beta=adapt_tau_2*v_t_beta+(1-adapt_tau_2)*(grad_LB_beta.^2);
    mt_hat_beta=m_t_beta./(1-(adapt_tau_1.^i));
    vt_hat_beta=v_t_beta./(1-(adapt_tau_2.^i));
    
    lambda_betaloading=[mu_VB_beta;func_v(B_VB_beta);d_VB_beta;gam_YJ_beta_transform];
    lambda_betaloading=lambda_betaloading+adapt_alpha*(mt_hat_beta./(sqrt(vt_hat_beta)+adapt_epsilon));
    mu_VB_beta=lambda_betaloading(1:length_mu_VB_beta,1);
    temp_B_VB_beta=lambda_betaloading(length_mu_VB_beta+1:length_mu_VB_beta+length_B_VB_beta,1);
    d_VB_beta=lambda_betaloading(length_mu_VB_beta+length_B_VB_beta+1:length_mu_VB_beta+length_B_VB_beta+length_d_VB_beta,1);
    gam_YJ_beta_transform=lambda_betaloading(length_mu_VB_beta+length_B_VB_beta+length_d_VB_beta+1:length_mu_VB_beta+length_B_VB_beta+length_d_VB_beta+length_gam_YJ_beta_transform,1);
    
    for j=1:num_factor_VB     
        B_VB_beta(:,j)=[zeros(j-1,1);temp_B_VB_beta(((j-1)*num_param_betaloading-sum(0:(j-2)))+1:j*num_param_betaloading-sum(0:(j-1)),1)];   
    end
    
    end
    end
    
    
    
    
end
%prediction

