function [indx]=obtain_index_C1(num_param)

count=1;
lag=0;
for i=1:num_param
    for j=i-lag:i
        indx(count,1)=i;
        indx(count,2)=j;
        count=count+1;
    end
    lag=lag+1;
end
indx=indx(1:end,:);
end