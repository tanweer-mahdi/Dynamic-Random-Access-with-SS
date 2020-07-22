function refined = stat_refine(data,percent)
%% This function weeds out outliers and NaN from output data matrix to provide meaningful results
n = size(data,2);
k = round(n*percent);
refined = data;
for i = 1:size(data,1)
    tmp = refined(i,:);
    [val,ind] = maxk(tmp,k);
    refined(:,ind) = [] ; %removing the outliers
    tmp = refined(i,:); %updating temp
    refined(:,find(isnan(tmp)==1)) = []; %removing the outliers
end
