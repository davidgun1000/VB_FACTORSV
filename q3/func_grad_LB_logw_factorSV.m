function [grad_LB_idio,grad_LB_fact,grad_LB_beta,grad_LB_factscores,log_weight]=func_grad_LB_logw_factorSV(y,mu_G_idio,C_G_idio,d_L_idio,D_L_idio,f_L_idio,F_L_idio,num_G_idio,num_L_idio,...
    C_L_star_idio,indx_CL_idio,indx_diag_CL_idio,func_v_eye_numL_idio,repmat_01_idio,...
    mu_G_fact,C_G_fact,d_L_fact,D_L_fact,f_L_fact,F_L_fact,num_G_fact,num_L_fact,C_L_star_fact,indx_CL_fact,indx_diag_CL_fact,func_v_eye_numL_fact,repmat_01_fact,...
    mu_VB_factscores,d_VB_factscores,gam_YJ_VB_factscores_transform,mu_VB_beta,B_VB_beta,d_VB_beta,gam_YJ_beta_transform,prior,dim_y,num_fact,num_param_factscores,num_param_betaloading,num_factor_VB)

    for j=1:num_fact
        epsilon_factscores{j,1}=randn(num_param_factscores,1);
        psi_theta_factscores{j,1}=mu_VB_factscores{j,1}+d_VB_factscores{j,1}.*epsilon_factscores{j,1};
        gam_YJ_VB_factscores{j,1}=logitcdf_gamYJ(gam_YJ_VB_factscores_transform{j,1});
        theta_factscores{j,1}=YJ_psi_to_theta(psi_theta_factscores{j,1},gam_YJ_VB_factscores{j,1});
    end

    s_VB_beta=randn(num_factor_VB,1);
    epsilon_VB_beta=randn(num_param_betaloading,1);
    psi_beta=mu_VB_beta+B_VB_beta*s_VB_beta+d_VB_beta.*epsilon_VB_beta;
    gam_YJ_beta=logitcdf_gamYJ(gam_YJ_beta_transform);
    theta_beta=YJ_psi_to_theta(psi_beta,gam_YJ_beta);
    
    for j=1:dim_y
        s_idio{j,1}=randn(num_L_idio+num_G_idio,1);
        s_G_idio{j,1}=s_idio{j,1}(1:num_G_idio,1);
        s_L_idio{j,1}=s_idio{j,1}(num_G_idio+1:end,1);
        var1_idio=((C_G_idio{j,1}')\s_G_idio{j,1});
        theta_G_idio{j,1}=mu_G_idio{j,1}+var1_idio;
        
        v_C_L_star_idio=f_L_idio{j,1}+F_L_idio{j,1}*theta_G_idio{j,1};
        id=v_C_L_star_idio~=0;
        v_C_L_temp_idio=v_C_L_star_idio(id,1);
        
        C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_CL_idio(:,1)'],[indx_CL_idio(:,2)']))=v_C_L_temp_idio;
        C_L_idio{j,1}=C_L_star_idio{j,1};
        temp_C_L_star_idio=exp(C_L_star_idio{j,1}(sub2ind(size(C_L_star_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)'])));
        C_L_idio{j,1}(sub2ind(size(C_L_idio{j,1}),[indx_diag_CL_idio(:,1)'],[indx_diag_CL_idio(:,2)']))=temp_C_L_star_idio';
        
        var2_idio=((C_L_idio{j,1}')\s_L_idio{j,1});
        var3_idio=((C_L_idio{j,1}')\D_L_idio{j,1});
        var3_var1_idio=var3_idio*var1_idio;
        mu_L_idio=d_L_idio{j,1}-var3_var1_idio;
        
        theta_L_idio{j,1}=mu_L_idio+var2_idio;
        
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
        temp_temp_idio=(var2_idio(indx_CL_idio(:,1),1)-var3_var1_idio(indx_CL_idio(:,1),1)).*(temp_GL_CL_idio(1,indx_CL_idio(:,2)))';
        
        grad_thetaL_G2_idio=-D2_star_idio.*temp_temp_idio;
        grad_mu_thetaL_G2_idio=F_L_idio{j,1}'*grad_thetaL_G2_idio;
        
        grad_LB_temp1_idio=G_G_idio+grad_mu_thetaL_G2_idio;
        grad_LB_temp2_idio=-D1_star_idio.*func_v(var1_idio*(((G_G_idio+grad_mu_thetaL_G2_idio-(D_L_idio{j,1}')*var4_idio)')/(C_G_idio{j,1}')));
        grad_LB_temp3_idio=G_L_idio;
        grad_LB_temp4_idio=-vec(var4_idio*var1_idio');
        grad_LB_temp5_idio=grad_thetaL_G2_idio;
        grad_LB_temp6_idio=vec(grad_thetaL_G2_idio*(theta_G_idio{j,1}'));
        
        grad_LB_idio(:,j)=[grad_LB_temp1_idio;grad_LB_temp2_idio;grad_LB_temp3_idio;grad_LB_temp4_idio;grad_LB_temp5_idio;grad_LB_temp6_idio];
    end
   
    for j=1:num_fact
        s_fact{j,1}=randn(num_L_fact+num_G_fact,1);
        s_G_fact{j,1}=s_fact{j,1}(1:num_G_fact,1);
        s_L_fact{j,1}=s_fact{j,1}(num_G_fact+1:end,1);
        var1_fact=((C_G_fact{j,1}')\s_G_fact{j,1});
        theta_G_fact{j,1}=mu_G_fact{j,1}+var1_fact;
        
        v_C_L_star_fact=f_L_fact{j,1}+F_L_fact{j,1}*theta_G_fact{j,1};
        id=v_C_L_star_fact~=0;
        v_C_L_temp_fact=v_C_L_star_fact(id,1);
        
        C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_CL_fact(:,1)'],[indx_CL_fact(:,2)']))=v_C_L_temp_fact;
        C_L_fact{j,1}=C_L_star_fact{j,1};
        temp_C_L_star_fact=exp(C_L_star_fact{j,1}(sub2ind(size(C_L_star_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)'])));
        C_L_fact{j,1}(sub2ind(size(C_L_fact{j,1}),[indx_diag_CL_fact(:,1)'],[indx_diag_CL_fact(:,2)']))=temp_C_L_star_fact';
        
        var2_fact=((C_L_fact{j,1}')\s_L_fact{j,1});
        var3_fact=((C_L_fact{j,1}')\D_L_fact{j,1});
        var3_var1_fact=var3_fact*var1_fact;
        mu_L_fact=d_L_fact{j,1}-var3_var1_fact;
        
        theta_L_fact{j,1}=mu_L_fact+var2_fact;
        D1_star_fact=func_v(func_dg(C_G_fact{j,1})+ones(num_G_fact,num_G_fact)-speye(num_G_fact));
        
        diag_CL_fact=diag(C_L_fact{j,1});
        diag_CL_kron_fact=kron(diag_CL_fact,[1;0]);
        diag_CL_kron_fact=diag_CL_kron_fact(1:2*num_L_fact-1,1);
        
        D2_star_fact=diag_CL_kron_fact+repmat_01_fact;
        [grad_param_states_fact]=obtain_grad_param_states_prior_fact(theta_G_fact{j,1},theta_L_fact{j,1},y,theta_factscores,theta_beta,num_G_fact,num_L_fact,j,prior,dim_y,num_fact);
        [grad_logq_fact]=obtain_grad_logq_SV(d_L_fact{j,1},D_L_fact{j,1},F_L_fact{j,1},C_L_fact{j,1},C_G_fact{j,1},theta_L_fact{j,1},s_G_fact{j,1},s_L_fact{j,1},...
            D2_star_fact,indx_CL_fact,func_v_eye_numL_fact);
        G_fact=grad_param_states_fact-grad_logq_fact;
        G_G_fact=G_fact(1:num_G_fact,1);
        G_L_fact=G_fact(num_G_fact+1:end,1);
        
        var4_fact=((C_L_fact{j,1})\G_L_fact);
        temp_GL_CL_fact=((G_L_fact')/(C_L_fact{j,1}'));
        temp_temp_fact=(var2_fact(indx_CL_fact(:,1),1)-var3_var1_fact(indx_CL_fact(:,1),1)).*(temp_GL_CL_fact(1,indx_CL_fact(:,2)))';
        
        grad_thetaL_G2_fact=-D2_star_fact.*temp_temp_fact;
        grad_mu_thetaL_G2_fact=F_L_fact{j,1}'*grad_thetaL_G2_fact;
        
        grad_LB_temp1_fact=G_G_fact+grad_mu_thetaL_G2_fact;
        grad_LB_temp2_fact=-D1_star_fact.*func_v(var1_fact*(((G_G_fact+grad_mu_thetaL_G2_fact-(D_L_fact{j,1}')*var4_fact)')/(C_G_fact{j,1}')));
        grad_LB_temp3_fact=G_L_fact;
        grad_LB_temp4_fact=-vec(var4_fact*var1_fact');
        grad_LB_temp5_fact=grad_thetaL_G2_fact;
        grad_LB_temp6_fact=vec(grad_thetaL_G2_fact*(theta_G_fact{j,1}'));
        
        grad_LB_fact(:,j)=[grad_LB_temp1_fact;grad_LB_temp2_fact;grad_LB_temp3_fact;grad_LB_temp4_fact;grad_LB_temp5_fact;grad_LB_temp6_fact];
    end
    
    ctraj_idio=[];
    ctraj_fact=[];
    for j=1:dim_y
        ctraj_idio=[ctraj_idio;theta_L_idio{j,1}'];
    end
    
    for j=1:num_fact
        ctraj_fact=[ctraj_fact;theta_L_fact{j,1}'];
    end
    
    [grad_betaloading,grad_factscores]=obtain_grad_factor_beta_loading_SV(theta_G_idio,theta_G_fact,theta_beta,theta_factscores,y,...
        ctraj_idio,ctraj_fact,num_param_factscores,num_param_betaloading,dim_y,num_fact);
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
    
    for j=1:num_fact
        [grad_logq_factscores]=obtain_grad_logq_factscores(theta_factscores{j,1},psi_theta_factscores{j,1},mu_VB_factscores{j,1},d_VB_factscores{j,1},gam_YJ_VB_factscores{j,1});
         diff_gradlogpos_gradlogq_factscores=grad_factscores(j,:)'-grad_logq_factscores;
         temp_deriv_theta_to_psi_factscores=deriv_theta_to_psi(psi_theta_factscores{j,1},gam_YJ_VB_factscores{j,1});
         grad_mu_VB_factscores=temp_deriv_theta_to_psi_factscores.*diff_gradlogpos_gradlogq_factscores;
         grad_d_VB_factscores=(temp_deriv_theta_to_psi_factscores.*epsilon_factscores{j,1}).*diff_gradlogpos_gradlogq_factscores;
         grad_gamYJ_factscores_transform=(deriv_theta_to_gamYJ(theta_factscores{j,1},psi_theta_factscores{j,1},gam_YJ_VB_factscores{j,1})).*diff_gradlogpos_gradlogq_factscores.*grad_gamYJ_to_gamYJtransform(gam_YJ_VB_factscores_transform{j,1});
         grad_LB_factscores(:,j)=[grad_mu_VB_factscores;grad_d_VB_factscores;grad_gamYJ_factscores_transform]; 
        
    end
    
    [log_weight]=compute_IW_factorSV(y,theta_G_idio,theta_L_idio,mu_G_idio,C_G_idio,mu_L_idio,C_L_idio,num_G_idio,num_L_idio,s_G_idio,s_L_idio,...
        theta_G_fact,theta_L_fact,mu_G_fact,C_G_fact,mu_L_fact,C_L_fact,num_G_fact,num_L_fact,s_G_fact,s_L_fact,...
        theta_beta,psi_beta,s_VB_beta,epsilon_VB_beta,mu_VB_beta,B_VB_beta,d_VB_beta,gam_YJ_beta,...
        theta_factscores,psi_theta_factscores,epsilon_factscores,mu_VB_factscores,d_VB_factscores,gam_YJ_VB_factscores,...
        dim_y,num_fact,num_param_factscores,num_param_betaloading,prior,num_factor_VB);
    
end