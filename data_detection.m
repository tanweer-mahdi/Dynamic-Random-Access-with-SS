function data_matrix = data_detection(XX,act_mat,channels,constellation)
data_matrix = zeros(size(XX));
XX =XX.*act_mat;
for i = 1:size(XX,1)
    pot = channels(i)*constellation;
    tmp = find(XX(i,:)~=0);
    data_matrix(i,tmp(1)) = 1; %pilot
    for ii= 2:length(tmp)
        [~,symbol] = min(abs(pot-XX(i,tmp(ii))));
        data_matrix(i,tmp(ii)) = constellation(symbol);
    end
end
end