function auset = spice_aud(snaps,phi,lambda)

%%
% This function detects users in random access opportunity using SPARSE
% ITERATIVE COVARIANCE ESTIMATION
% snaps = data snapshots
% phi = spreading sequences
% lambda = threshold for user decision
% auset = estimated set of active users

%%

M = size(snaps,1);
N = size(phi,2);
Rhat = snaps*snaps'/size(snaps,2);
aug_phi = [phi(1:M,:) eye(M)];
s_hat     = fft(snaps,N,1);
pp         = sum(s_hat.*conj(s_hat),2)/(size(snaps,2)*M^2);
p = [pp; diag(Rhat)];
sigma = trace(Rhat)/M; % initial estimate of noise variance
% initializing R
ff = fft(p(1:N));
ff(1) = ff(1) + sigma;
R = toeplitz((ff(1:M)));
w = vecnorm(aug_phi).^2/trace(Rhat);
w = abs(diag(aug_phi'*inv(Rhat)*aug_phi))'/M;
gamma  = sum(w(N+1:N+M));
Rhat = sqrtm(Rhat);
for iter=1:5
    temp1 = R\Rhat;
    temp2 = norm(temp1,'fro');
    temp3 = vecnorm((phi(1:M,:)'*temp1)'); % row vector
    rho = temp3*(p(1:N).*sqrt(w(1:N)')) + sqrt(gamma)*sigma*temp2;
    sigma = sigma*temp2/(sqrt(gamma)*rho);
    p(1:N) = (p(1:N)./(sqrt(w(1:N)')*rho)).*temp3';
    p(N+1:end) = sigma;
    %update R
    ff = fft(p(1:N));
    ff(1) = ff(1) + sigma;
    R = toeplitz((ff(1:M)));
end
auset = sort(find(log(abs(p(1:N)))>lambda))';
end
