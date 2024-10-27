%parpool(28)
load('DATA_SP100.mat');
load('InitValuesFactor.mat');
T_est=1000;
y=y_used_demean;
dim_y=size(y,1);
num_fact=4;
num_MonteCarlo=1000;
psi=2.75*ones(dim_y+num_fact,1);
sig=0.2*ones(dim_y+num_fact,1);
kapha=-0.6*ones(dim_y,1);

phi=exp(psi)./(1+exp(psi));
alpha=log(exp(sig)-1);
step_ahead=10;

for i=1:num_fact
    B(:,i)=[zeros(i-1,1);1;zeros(dim_y-i,1)];
end
burn=1000;
nit=40000;
iter=burn+nit;
N=100;
T=length(y(1,:)');

%y_predict=y(:,T+1);
%y=y(:,1:T);
D1=3;
for i=1:dim_y
    scale_idio(i,1)=1;
    V1(:,:,i)=0.01*eye(D1);
end
target_accept=0.20;
D2=2;
for i=1:num_fact
    scale_factor(i,1)=1;
    V2(:,:,i)=0.01*eye(D2);
end
%fact_score=randn(num_fact,T);
fact_score=[mu_VB_factscores{1,1}';mu_VB_factscores{2,1}';mu_VB_factscores{3,1}';mu_VB_factscores{4,1}';];
%ctraj=randn(dim_y+num_fact,T);
 parfor s=1:dim_y
   [X,W,A,~]=smc_mu_factorSV_noncentered(y(s,:),B(s,:),fact_score,phi(s,1),sig(s,1),kapha(s,1),N,T,s,dim_y);
   ctraj(s,:)=ancestor_trace(X,W,A);
 end
 parfor s=dim_y+1:dim_y+num_fact
   [X,W,A,~]=smc_mu_factorSV_noncentered(0,0,fact_score,phi(s,1),sig(s,1),0,N,T,s,dim_y);
   ctraj(s,:)=ancestor_trace(X,W,A);
 end

for t=1:T
        ctraj_idio_temp(:,t)=ctraj(1:dim_y,t);
        ctraj_fact_temp(:,t)=ctraj(dim_y+1:dim_y+num_fact,t);
        Vt_inv=diag(1./(exp(sig(1:dim_y,1).*ctraj_idio_temp(:,t)+kapha(1:dim_y,1))));
        Dt_inv=diag(exp(-(sig(dim_y+1:dim_y+num_fact,1).*ctraj_fact_temp(:,t))));
        
        %var_ft=inv((B'*Vt_inv*B)+Dt_inv);
        %[var_ft]=jitChol(var_ft);
        
        [chol_var_ft,flag]=chol(inv((B'*Vt_inv*B)+Dt_inv),'lower');
        if flag==0
        %chol_var_ft=chol(var_ft);
        mean_ft=(chol_var_ft*chol_var_ft')*B'*(Vt_inv*y(:,t));
        fact_score(:,t)=mvnrnd(mean_ft,chol_var_ft*chol_var_ft');      
    
        end
end 
 
 
 
B0_deep=1;
hp_sig2=10;
hp_B=1;
a0_phi=20;
b0_phi=1.5;

accept1=zeros(num_fact,1);
accept2=zeros(dim_y+num_fact,1);
hp_mu=4;

for i=1:dim_y
    ctraj_idio_store{i,1}=zeros(1,T);
end

for i=1:num_fact
    ctraj_factor_store{i,1}=zeros(1,T);
    latent_factor_store{i,1}=zeros(1,T);
end


w_portfolio=(1/dim_y)*ones(dim_y,1);
count=0;
for i=1:iter
        i
        %tic
%         [phi(91,1),phi(92,1),phi(93,1),phi(94,1)]
%         [B(1,1),B(2,2),B(3,3),B(4,4)]
%         [ctraj(dim_y+num_fact,1:100:1000)]
        
    for s=1:dim_y
        s_bar=min(s,num_fact);
        F=fact_score(1:s_bar,:)';
        arg_Vi_inv=1./(exp(sig(s,1).*ctraj(s,:)+kapha(s,1)));
        Vi_inv=diag(arg_Vi_inv);
        B_var=inv((F'*Vi_inv*F)+eye(s_bar));
        [B_var]=jitChol(B_var);
        chol_B=chol(B_var);
        B_mean=B_var*F'*((Vi_inv*y(s,:)'));
        temp_B=(mvnrnd(B_mean',chol_B'*chol_B));        
    
        
        
            
        if sum(isnan(temp_B))>0 | sum(isinf(temp_B))>0
           B(s,1:s_bar)=B(s,1:s_bar); 
            
        else
           if s<=num_fact
                
               if temp_B(1,s)<=0
                  B(s,1:s_bar)=B(s,1:s_bar);
               else
                  B(s,1:s_bar)=temp_B; 
               end
           else
              B(s,1:s_bar)=temp_B; 
           end 
        end           
        
    end
      
      for s=1:num_fact
          ctraj_fact_temp2(s,:)=sig(dim_y+s,1).*ctraj(dim_y+s,:);  
      end
      
      
      k1=dim_y;
      for s=1:num_fact
        ind=find(abs(B(:,s))==max(abs(B(:,s))));
        %max_B=max(B(s:dim_y,s));
        max_B=B(s,s);
        %max_B=abs(B(s,s));
        B_star=B(s:dim_y,s)./max_B;
        fact_star=max_B.*fact_score(s,:);
        
        ctraj_trans=ctraj_fact_temp2(s,:)+2*log(abs(max_B));
        mu_fact_old=log(max_B^2);
        mean_prop=(sum(ctraj_trans(1,2:T-1))+((ctraj_trans(1,T)-phi(dim_y+s,1)*ctraj_trans(1,1))./(1-phi(dim_y+s,1))))./(T-1+1/B0_deep);
        var_prop=((sig(dim_y+s,1)^2)/((1-phi(dim_y+s,1))^2))/(T-1+1/B0_deep);
        mu_fact_new=normrnd(mean_prop,sqrt(var_prop));
        if sum(isinf(mu_fact_new))>0 | sum(isnan(mu_fact_new))>0 | sum(isnan(B_star))>0 | sum(isinf(B_star))>0
        else
        A2=rand();
        [chol_comp1_old,flag]=chol(hp_B*exp(-mu_fact_old)*eye(k1),'lower');
        if flag==0
        comp1_old=logmvnpdf(B_star',zeros(1,k1),chol_comp1_old*chol_comp1_old');
        comp2_old=logmvnpdf(ctraj_trans(1,1),mu_fact_old,(sig(dim_y+s,1)^2)/(1-phi(dim_y+s,1)^2));       
        comp3_old=log(exp(mu_fact_old/2-exp(mu_fact_old)/(2*hp_B)));
        comp4_old=logmvnpdf(mu_fact_old,0,B0_deep*(sig(dim_y+s,1)^2)/((1-phi(dim_y+s,1))^2));
        [chol_comp1_new,flag]=chol(hp_B*exp(-mu_fact_new)*eye(k1),'lower');
        if flag==0
        comp1_new=logmvnpdf(B_star',zeros(1,k1),chol_comp1_new*chol_comp1_new');
        comp2_new=logmvnpdf(ctraj_trans(1,1),mu_fact_new,(sig(dim_y+s,1)^2)/(1-phi(dim_y+s,1)^2));
        comp3_new=log(exp(mu_fact_new/2-exp(mu_fact_new)/(2*hp_B)));
        comp4_new=logmvnpdf(mu_fact_new,0,B0_deep*(sig(dim_y+s,1)^2)/((1-phi(dim_y+s,1))^2));
        R2=exp(comp1_new+comp2_new+comp3_new-comp1_old-comp2_old-comp3_old+comp4_old-comp4_new);
        C2=min(1,R2);
        if A2<=C2
           lam1=exp(mu_fact_new/2);
           accept1(s,1)=accept1(s,1)+1;   
        else
           lam1=max_B;
        end     
        B(:,s)=(lam1/max_B).*B(:,s);
        fact_score(s,:)=(max_B/lam1).*fact_score(s,:);
        ctraj_fact_temp2(s,:)=ctraj_fact_temp2(s,:)+2*log(abs(max_B/lam1));
        
        end
        
        end
        end
        k1=k1-1;
      end
      
    for s=1:num_fact  
    ctraj(dim_y+s,:)=ctraj_fact_temp2(s,:)./sig(dim_y+s,1);
    end
%     

    %sampling the latent factor f
    for t=1:T
        ctraj_idio_temp(:,t)=ctraj(1:dim_y,t);
        ctraj_fact_temp(:,t)=ctraj(dim_y+1:dim_y+num_fact,t);
        Vt_inv=diag(1./(exp(sig(1:dim_y,1).*ctraj_idio_temp(:,t)+kapha(1:dim_y,1))));
        Dt_inv=diag(exp(-(sig(dim_y+1:dim_y+num_fact,1).*ctraj_fact_temp(:,t))));
        
        %var_ft=inv((B'*Vt_inv*B)+Dt_inv);
        %[var_ft]=jitChol(var_ft);
        
        [chol_var_ft,flag]=chol(inv((B'*Vt_inv*B)+Dt_inv),'lower');
        if flag==0
        %chol_var_ft=chol(var_ft);
        mean_ft=(chol_var_ft*chol_var_ft')*B'*(Vt_inv*y(:,t));
        fact_score(:,t)=mvnrnd(mean_ft,chol_var_ft*chol_var_ft');      
    
        end
    end
    
    parfor s=1:dim_y
        [u1]=obtain_random_numbers(phi(s,1),sig(s,1),kapha(s,1),ctraj(s,:),N,T); 
        [X,W,A,lik,u1_res_rand]=csmc_mu_factorSV_corr_noncentered(y(s,:),B(s,:),fact_score,phi(s,1),sig(s,1),kapha(s,1),N,T,ctraj(s,:),u1,s,dim_y); 
        temp_param=[kapha(s,1),psi(s,1),alpha(s,1)];
        R1=mvnrnd(temp_param,scale_idio(s,1).*V1);
        kapha_star=R1(1,1);
        psi_star=R1(1,2);
        alpha_star=R1(1,3);
        phi_star=exp(psi_star)/(1+exp(psi_star));
        sig_star=log(exp(alpha_star)+1);
        
        [X_star,W_star,A_star,lik_star]=smc_mu_factorSV_corr_noncentered(y(s,:),B(s,:),fact_score,phi_star,sig_star,kapha_star,N,T,u1,u1_res_rand,s,dim_y);
        
        if sum(sum(isnan(W)))>0 | sum(sum(isnan(W_star)))>0 | sum(sum(isnan(X)))>0 | sum(sum(isnan(X_star)))>0 | phi_star>=0.999
           psi(s,1)=psi(s,1);
           alpha(s,1)=alpha(s,1);
           kapha(s,1)=kapha(s,1);
           phi(s,1)=phi(s,1);
           sig(s,1)=sig(s,1);
           ctraj(s,:)=ctraj(s,:);
           thetasave(i,:,s)=[kapha(s,1),psi(s,1),alpha(s,1)];
        else
        
        prior_kapha=log(normpdf(kapha(s,1),0,sqrt(hp_sig2)));
        prior_psi=(a0_phi-1).*log((1+phi(s,1))./2)+(b0_phi-1).*log((1-phi(s,1))./2)+((exp(psi(s,1)))./((1+exp(psi(s,1))).^2));
        prior_alpha=log(2)-log(pi)-log(1+sig(s,1).^2)+(exp(alpha(s,1))./(1+exp(alpha(s,1))));
        prior=prior_kapha+prior_psi+prior_alpha;
        post=prior+lik;
        
        prior_kapha_star=log(normpdf(kapha_star,0,sqrt(hp_sig2)));
        prior_psi_star=(a0_phi-1).*log((1+phi_star)./2)+(b0_phi-1).*log((1-phi_star)./2)+((exp(psi_star))./((1+exp(psi_star)).^2));
        prior_alpha_star=log(2)-log(pi)-log(1+sig_star.^2)+(exp(alpha_star)./(1+exp(alpha_star)));
        prior_star=prior_kapha_star+prior_psi_star+prior_alpha_star;
        post_star=lik_star+prior_star;
        
        r1 = exp(post_star-post);
        C1=min(1,r1);
        A1=rand();
        if A1<=C1
           psi(s,1)=psi_star;
           alpha(s,1)=alpha_star;
           kapha(s,1)=kapha_star;
           X=X_star;
           W=W_star;
           A=A_star;
           phi(s,1)=phi_star;
           sig(s,1)=sig_star; 
        end
        thetasave(i,:,s)=[kapha(s,1),psi(s,1),alpha(s,1)];
        ctraj(s,:)=backward_simulation_SV_noncentered(X,W,A,phi(s,1),sig(s,1),kapha(s,1));
        if i>100
           scale_idio(s,1)=update_sigma(scale_idio(s,1),C1,target_accept,i,3); 
        end   

        end
    end
    
    for s=1:dim_y
         if i>100
            V1(:,:,s)=cov(thetasave(1:i,:,s));
            V1(:,:,s)=jitChol(V1(:,:,s));
         end
    end
 
    parfor s=dim_y+1:dim_y+num_fact
       [u1]=obtain_random_numbers(phi(s,1),sig(s,1),0,ctraj(s,:),N,T);
       [X,W,A,lik,u1_res_rand]=csmc_mu_factorSV_corr_noncentered(0,0,fact_score,phi(s,1),sig(s,1),0,N,T,ctraj(s,:),u1,s,dim_y);
       temp_param=[psi(s,1),alpha(s,1)];
       R1=mvnrnd(temp_param,scale_factor(s-dim_y,1).*V2(:,:,s-dim_y));
       psi_star=R1(1,1);
       alpha_star=R1(1,2);
       phi_star=exp(psi_star)/(1+exp(psi_star));
       sig_star=log(exp(alpha_star)+1);
       
       [X_star,W_star,A_star,lik_star]=smc_mu_factorSV_corr_noncentered(0,0,fact_score,phi_star,sig_star,0,N,T,u1,u1_res_rand,s,dim_y);
       
       if sum(sum(isnan(W)))>0 | sum(sum(isnan(W_star)))>0 | sum(sum(isnan(X)))>0 | sum(sum(isnan(X_star)))>0 | phi_star>=0.995
          psi(s,1)=psi(s,1);
          alpha(s,1)=alpha(s,1);
          phi(s,1)=phi(s,1);
          sig(s,1)=sig(s,1);
          ctraj(s,:)=ctraj(s,:);
          thetasave2(i,:,s-dim_y)=[psi(s,1),alpha(s,1)];
       else
       prior_psi=(a0_phi-1).*log((1+phi(s,1))./2)+(b0_phi-1).*log((1-phi(s,1))./2)+((exp(psi(s,1)))./((1+exp(psi(s,1))).^2));
       prior_alpha=log(2)-log(pi)-log(1+sig(s,1).^2)+(exp(alpha(s,1))./(1+exp(alpha(s,1))));
       prior=prior_psi+prior_alpha;
       post=prior+lik;
       
       prior_psi_star=(a0_phi-1).*log((1+phi_star)./2)+(b0_phi-1).*log((1-phi_star)./2)+((exp(psi_star))./((1+exp(psi_star)).^2));
       prior_alpha_star=log(2)-log(pi)-log(1+sig_star.^2)+(exp(alpha_star)./(1+exp(alpha_star)));
       prior_star=prior_psi_star+prior_alpha_star;
       post_star=lik_star+prior_star;
       
       r1 = exp(post_star-post);
       C1=min(1,r1);
       A1=rand();
       if A1<=C1
          psi(s,1)=psi_star;
          alpha(s,1)=alpha_star;
          X=X_star;
          W=W_star;
          A=A_star;
          phi(s,1)=phi_star;
          sig(s,1)=sig_star; 
       end
       thetasave2(i,:,s-dim_y)=[psi(s,1),alpha(s,1)];
       ctraj(s,:)=backward_simulation_SV_noncentered(X,W,A,phi(s,1),sig(s,1),0); 
       if i>100
          scale_factor(s-dim_y,1)=update_sigma(scale_factor(s-dim_y,1),C1,target_accept,i,2); 
       end
       
       end
       
    end
    
    for s=dim_y+1:dim_y+num_fact
         if i>100
            V2(:,:,s-dim_y)=cov(thetasave2(1:i,:,s-dim_y));
            V2(:,:,s-dim_y)=jitChol(V2(:,:,s-dim_y)); 
         end 
    end
 
    if i>burn
       for j=1:dim_y
           ctraj_idio_store{j,1}=ctraj_idio_store{j,1}+ctraj(j,:);
       end
       
       for j=1:num_fact
           ctraj_factor_store{j,1}=ctraj_factor_store{j,1}+ctraj(dim_y+j,:);
           latent_factor_store{j,1}=latent_factor_store{j,1}+fact_score(j,:);
       end

    end
   

     

    %computing the predictive density
     
     if i>burn & mod(i,20)==0
%         ctraj_factTplusone=ctraj(dim_y+1:end,T);
%         ctraj_idioTplusone=ctraj(1:dim_y,T);
        count=count+1;
%         for k=1:step_ahead
%             ctraj_factTplusone=phi(dim_y+1:end,1).*ctraj_factTplusone+randn(num_fact,1);
%             ctraj_idioTplusone=phi(1:dim_y,1).*ctraj_idioTplusone+randn(dim_y,1);
%             biglambda = diag([exp(sig(dim_y+1).*ctraj_factTplusone)]);
%             bigV=diag([exp(sig(1:dim_y,1).*ctraj_idioTplusone+kapha(1:dim_y,1))]);
%             covmat=B*biglambda*B'+bigV;
%             y_pred(count,:,k)=(mvnrnd(zeros(1,dim_y),covmat))';
%             y_pred_portfolio(count,:,k)=sum(w_portfolio'.*y_pred(count,:,k),2);
%             rhomat_pred{k,1}(:,:,count)=diag(diag(covmat))^-0.5 * (covmat * diag(diag(covmat))^-0.5);
%         end
        
        for j=1:num_fact
            Post_param.B{j,1}(count,:)=B(j:end,j)';
        end
    %    Post_param.B1(count,:)=B(:,1)';
    %     Post_param.B2(i,:)=B(:,2)';
    %     Post_param.B3(i,:)=B(:,3)';
    %     Post_param.B4(i,:)=B(:,4)';

       Post_param.sig(count,:)=sig';
       Post_param.phi(count,:)=phi';
       Post_param.kapha(count,:)=kapha';
       Post_param.psi(count,:)=psi';
       Post_param.alpha(count,:)=alpha';

       id=(1:1:100)*10;
       
       for j=1:dim_y
           Post_latent.ctraj_idio{j,1}(count,:)=ctraj(j,id);
       end
       for j=1:num_fact
           Post_latent.ctraj_factor{j,1}(count,:)=ctraj(dim_y+j,id);
           Post_latent.latent_factor{j,1}(count,:)=fact_score(j,id);
       end
       
%        Post_latent.ctraj1(count,:)=ctraj(1,id);
%        Post_latent.ctraj2(count,:)=ctraj(2,id);
%        Post_latent.ctraj3(count,:)=ctraj(3,id);
%        Post_latent.ctraj4(count,:)=ctraj(4,id);
%        Post_latent.ctraj5(count,:)=ctraj(5,id);
%        Post_latent.ctraj6(count,:)=ctraj(6,id);
%        Post_latent.ctraj7(count,:)=ctraj(7,id);
%        Post_latent.ctraj8(count,:)=ctraj(8,id);
%        Post_latent.ctraj9(count,:)=ctraj(9,id);
%        Post_latent.ctraj10(count,:)=ctraj(10,id);
%        Post_latent.ctraj11(count,:)=ctraj(11,id);
%        Post_latent.ctraj12(count,:)=ctraj(12,id);
%        Post_latent.ctraj13(count,:)=ctraj(13,id);
%        Post_latent.ctraj14(count,:)=ctraj(14,id);
%        Post_latent.ctraj15(count,:)=ctraj(15,id);
%        Post_latent.ctraj16(count,:)=ctraj(16,id);
%        Post_latent.ctraj17(count,:)=ctraj(17,id);
%        Post_latent.ctraj18(count,:)=ctraj(18,id);
%        Post_latent.ctraj19(count,:)=ctraj(19,id);
%        Post_latent.ctraj20(count,:)=ctraj(20,id);
%        Post_latent.ctraj21(count,:)=ctraj(21,id);
%        Post_latent.ctraj22(count,:)=ctraj(22,id);
%        Post_latent.ctraj23(count,:)=ctraj(23,id);
%        Post_latent.ctraj24(count,:)=ctraj(24,id);
%        Post_latent.ctraj25(count,:)=ctraj(25,id);
%        Post_latent.ctraj26(count,:)=ctraj(26,id);
%         
%        Post_latent.ctrajf1(count,:)=ctraj(27,id);
%        Post_latent.fact_score1(count,:)=fact_score(1,id);
        
        
        
    end 
    if mod(i,2000)==0
       save('corr_SP100_noncentered_pred_4factor_deep.mat','Post_param','Post_latent');
       %save('corr_SP100_noncentered_pred_4factor.mat','Post_param','Post_latent','y_pred','y_pred_portfolio','rhomat_pred');
    
    end
    
    %toc
end

for j=1:dim_y
    Post.ctraj_idio_mean{j,1}=ctraj_idio_store{j,1}./nit; 
end
for j=1:num_fact
    Post.ctraj_factor_mean{j,1}=ctraj_factor_store{j,1}./nit;
    Post.latent_factor_mean{j,1}=latent_factor_store{j,1}./nit;
end
save('corr_SP100_noncentered_pred_4factor_deep.mat','Post_param','Post_latent','Post');
%save('corr_SP100_noncentered_pred_4factor.mat','Post_param','Post_latent','y_pred','y_pred_portfolio','rhomat_pred');

%     Post.ctraj1sum=ctraj1_store./nit;
%     Post.ctraj2sum=ctraj2_store./nit;
%     Post.ctraj3sum=ctraj3_store./nit;
%     Post.ctraj4sum=ctraj4_store./nit;
%     Post.ctraj5sum=ctraj5_store./nit;
%     Post.ctraj6sum=ctraj6_store./nit;
%     Post.ctraj7sum=ctraj7_store./nit;
%     Post.ctraj8sum=ctraj8_store./nit;
%     Post.ctraj9sum=ctraj9_store./nit;
%     Post.ctraj10sum=ctraj10_store./nit;
%     Post.ctraj11sum=ctraj11_store./nit;
%     Post.ctraj12sum=ctraj12_store./nit;
%     Post.ctraj13sum=ctraj13_store./nit;
%     Post.ctraj14sum=ctraj14_store./nit;
%     Post.ctraj15sum=ctraj15_store./nit;
%     Post.ctraj16sum=ctraj16_store./nit;
%     Post.ctraj17sum=ctraj17_store./nit;
%     Post.ctraj18sum=ctraj18_store./nit;
%     Post.ctraj19sum=ctraj19_store./nit;
%     Post.ctraj20sum=ctraj20_store./nit;
%     Post.ctraj21sum=ctraj21_store./nit;
%     Post.ctraj22sum=ctraj22_store./nit;
%     Post.ctraj23sum=ctraj23_store./nit;
%     Post.ctraj24sum=ctraj24_store./nit;
%     Post.ctraj25sum=ctraj25_store./nit;
%     Post.ctraj26sum=ctraj26_store./nit;
%     Post.ctrajf1sum=ctrajf1_store./nit;
% 
%     Post.fact1sim=fact1_store./nit;

