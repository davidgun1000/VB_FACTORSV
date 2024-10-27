function [indx]=obtain_index_C2(num_states)

count=1;
for i=1:num_states
    for j=i-1:i
        indx(count,1)=i;
        indx(count,2)=j;
        count=count+1;
    end
end
indx=indx(2:end,:);
% length_indx=length(indx);
% 
% lag=0;
% count=length(indx)+1;
% for i=num_states+1:num_states+num_param
%     for j=i-lag:i
%         indx(count,1)=i;
%         indx(count,2)=j;
%         count=count+1;
%     end
%     lag=lag+1;   
% end

end