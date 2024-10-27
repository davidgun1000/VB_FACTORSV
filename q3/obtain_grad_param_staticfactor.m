function [grad_beta,grad_factscores,grad_sig2]=obtain_grad_param_staticfactor(theta_beta,theta_factscores,theta_sig2,y,prior,dim_y,num_fact,num_param_beta,num_param_factscores,num_param_sig2)

    T=length(y(1,:)');
    sig2=exp(theta_sig2);
    beta_loading_temp=theta_beta;
    for i=1:num_fact
        fact_score(i,1:T)=theta_factscores{i,1};
        beta_loading(:,i)=[zeros(i-1,1);beta_loading_temp(((i-1)*dim_y-sum(0:(i-2)))+1:i*dim_y-sum(0:(i-1)),1)];
        beta_loading(i,i)=exp(beta_loading(i,i));
    end
        
    grad_beta_loading=zeros(dim_y,num_fact);
    temp_beta=zeros(dim_y,num_fact);
    for i=1:dim_y
        i_bar=min(i,num_fact);
        F=fact_score(1:i_bar,:)';
        grad_beta_loading(i,1:i_bar)=((y(i,:)-(F*beta_loading(i,1:i_bar)')').*(1./sig2(i,1)))*F-beta_loading(i,1:i_bar);
        if i<=num_fact
           temp_beta(i,1:i_bar)=((y(i,:)-(F*beta_loading(i,1:i_bar)')').*(1./sig2(i,1)))*F;
           save_beta=temp_beta(i,i)*beta_loading(i,i)+1-beta_loading(i,i)^2;
           grad_beta_loading(i,i)=save_beta;
        end
    end
    
    temp1=((y(:,1:T)-beta_loading*fact_score)./sig2);
    temp1_reshape=reshape(temp1,dim_y,1,T);
    temp1_reshape_transp=multitransp(temp1_reshape);
    
    beta_loading_repmat=repmat(beta_loading,1,1,T);
    arg1=multiprod(temp1_reshape_transp,beta_loading_repmat);
    arg1=reshape(arg1,num_fact,T);
    
    for k=1:num_fact
        arg2(k,:)=fact_score(k,:);
    end
    grad_fact_score=arg1-arg2;
    
    for i=1:dim_y
        grad_sig2_temp1(i,1)=(-(T/2)*(1/sig2(i,1))+0.5*(1./(sig2(i,1)^2))*sum((y(i,:)-beta_loading(i,:)*fact_score).^2))*sig2(i,1);  
        grad_sig2_temp2(i,1)=((-((prior.v0)/2)-1)*(1/sig2(i,1))+(prior.s0/2).*(1/(sig2(i,1)^2)))*sig2(i,1)+1;
        grad_sig2(i,1)=grad_sig2_temp1(i,1)+grad_sig2_temp2(i,1);
    end
    
    grad_beta=[];
    for i=1:num_fact
        grad_beta=[grad_beta;grad_beta_loading(i:end,i)];
        grad_factscores(i,:)=grad_fact_score(i,:);
    end
    
%     
%     grad_param=[];
%     for i=1:num_fact
%         grad_param=[grad_param;grad_beta_loading(i:end,i)];
%     end
%     
%     for i=1:num_fact
%         grad_param=[grad_param;grad_fact_score(i,:)'];
%     end
%     
%     grad_param=[grad_param;grad_sig2];
    
    
end