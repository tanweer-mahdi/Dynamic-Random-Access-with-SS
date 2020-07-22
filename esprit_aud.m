function [auset varest] = esprit_aud(snaps,N)
m = size(snaps,1);
[U,E,V] = svd(snaps);
l = diag(E).^2/size(snaps,2); % Eigenvalues of Covariance matrix of snapshots
%% Bayesian Information Criteria
bic = zeros(1,m-1);
n = size(snaps,2);
for k=1:m-1
    bic(k) = -n*sum(log(l(k+1:m))) + n*(m-k)*log(mean(l(k+1:m))) + 0.5*k*(2*m-k)*log(n);
end
[~,uinit] = min(bic);
varest = mean(l(uinit+1:end))/2; %estimated noise variance
%% ESPRIT
omega = (2*pi)/N*(0:N-1);
U1 = U(1:end-1,1:uinit);
U2 = U(2:end,1:uinit);
X = U1\U2;
e = eig(X);
om = wrapTo2Pi(angle(e));
auset = [];

for i=1:length(om)
    [~,temp] = min(abs(om(i)-omega));
    auset = [auset temp];
end
auset = sort(unique(auset));
end