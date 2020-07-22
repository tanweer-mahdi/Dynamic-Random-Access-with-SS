function [y,rel] = channel_estimator(XX,act_mat,l)

% l = if l==1, then the function for dynamic access is chosen
% y = estimated channels of reliable users
% rel = list of reliable users
if l==1
nz = sum(act_mat,2); % number of non-zero entries in each row
%% Deciding which users are reliable 
XX = log(XX).*act_mat;
rel_ind = [];
rel_channel_gain = [];
for i=1:size(XX,1)
    tmp = real(XX(i,:));
    hh = exp(mean(tmp(tmp~=0)));
    eta = (hh/sqrt(nz(i)))*sqrt(norm(exp(tmp(tmp~=0)-log(hh)))^2 - 2*sum(exp(tmp(tmp~=0)-log(hh))) + nz(i));
    %eta = std(tmp(tmp~=0))
    if 2*hh*sin(pi/4) > 6*eta
        rel_ind = [rel_ind i];
        rel_channel_gain = [rel_channel_gain hh];
    end
end
%plot(real(XX)')
rel = rel_ind;
%% Estimating channel angles
% computing the candidate angles
cand = (sum(mod(4*imag(XX),2*pi),2)./nz)';
pot = [cand/4; cand/4 + 2*pi/4; cand/4 + 4*pi/4; cand/4 + 6*pi/4]; % potential angles
pilot = [];
for i=1:size(XX,1)
    tmp = find(abs(XX(i,:))>0);
    if isempty(tmp)
        pilot = [pilot 0]; % inserting zeros for all zero activation case
    else
    pilot = [pilot mod(imag(XX(i,tmp(1))),2*pi)];
    end
end
%pilot = mod(imag(XX(:,1)),2*pi);
rel_channel_angle = [];
for i=1:length(rel_ind)
    [~,ang] = min(abs(pot(:,rel_ind(i))-pilot(rel_ind(i))));
    rel_channel_angle = [rel_channel_angle pot(ang,rel_ind(i))];
end
y = rel_channel_gain.*exp(1i*rel_channel_angle); %estimated complex channel gains
end  
end