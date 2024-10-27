function [grad_factscores]=obtain_grad_factor_scores(theta_G_idio,theta_G_fact,theta_beta,theta_factscores,y,ctraj_idio,ctraj_fact,num_param_factscores,num_param_betaloading,dim_y,num_fact)

      T=length(y(1,:)');
      beta_loading_temp=theta_beta;
      for i=1:num_fact
          fact_score(i,1:T)=theta_factscores{i,1};
          beta_loading(:,i)=[zeros(i-1,1);beta_loading_temp(((i-1)*dim_y-sum(0:(i-2)))+1:i*dim_y-sum(0:(i-1)),1)];
          beta_loading(i,i)=exp(beta_loading(i,i));
      end
       
      for i=1:dim_y
          temp_kapha_idio(i,1)=theta_G_idio{i,1}(1,1);
          temp_sig_idio(i,1)=log(exp(theta_G_idio{i,1}(3,1))+1);
      end
      
      temp1=y(:,1:T)-beta_loading*fact_score;
      temp1_reshape=reshape(temp1,dim_y,1,T);
      temp1_reshape_transp=multitransp(temp1_reshape);
         
      ctraj_idio_inv_reshape=reshape(exp(-(temp_sig_idio.*ctraj_idio+temp_kapha_idio)),dim_y,1,T);
      ctraj_idio_inv_reshape_transp=multitransp(ctraj_idio_inv_reshape);
       
      temp2=temp1_reshape_transp.*ctraj_idio_inv_reshape_transp;
      beta_loading_repmat=repmat(beta_loading,1,1,T);
      arg1=multiprod(temp2,beta_loading_repmat);
      arg1=reshape(arg1,num_fact,T);
        
        
      for k=1:num_fact
          temp_sig_fact(k,1)=log(exp(theta_G_fact{k,1}(2,1))+1);
          arg2(k,:)=fact_score(k,:).*exp(-(temp_sig_fact(k,1).*ctraj_fact(k,1:T)));
      end
       
      grad_factscores_temp=arg1-arg2;
      
      for i=1:num_fact
          grad_factscores(i,:)=grad_factscores_temp(i,:);
      end
                 
end


