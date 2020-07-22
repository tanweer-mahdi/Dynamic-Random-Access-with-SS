function [auset_ref,act_mat,XX] = act_detect(auset,phi,yn,pfa,varest)
%% act_detect stands for Activity Detection
% auset = estimated set of users
% phi = set of spreading sequence
% yn = received signal
% act_mat = binary activity matrix
% XX = least square solution
phid = phi(:,auset);
X = inv(phid'*phid)*phid';
XX = X*yn;
act_mat = [];
for i=1:size(XX,1)
    suffstat = abs(XX(i,:)).^2;
    qn = X(i,:);
    % Calculation of ROC
    %             pfa = linspace(0.00001,0.1,1000);
    %             pd = zeros(1,length(pfa));
    %             for kk=1:length(pfa)
    thresh = gaminv(1-pfa,1,2*qn*qn'*varest);
    %             pd(kk) = 1-gamcdf(thresh,1,2*qn*qn'*varest + rp(auset(i)));
    %             end
    %             plot(pfa,pd)
    %             hold on
    %thresh = gaminv(1-pd,1,2*qn*qn'*varest + rp(auset(i))/2);
    act_mat = [act_mat; suffstat>thresh];
end
%% Refining estimated user set
nz = sum(act_mat,2);
if sum(nz>1)< length(auset)
    auset_ref = auset(find(nz>1));
    XX = phi(:,auset_ref)\yn;
    %act_mat = act_mat(find(nz>1),:);
    % Refining activation matrix
    phid = phi(:,auset_ref);
    X = inv(phid'*phid)*phid';
    XX = X*yn;
    act_mat = [];
    
    for i=1:size(XX,1)
        suffstat = abs(XX(i,:)).^2;
        qn = X(i,:);
        thresh = gaminv(1-pfa,1,2*qn*qn'*varest);
        act_mat = [act_mat; suffstat>thresh];
    end
else
    auset_ref = auset;
end

end