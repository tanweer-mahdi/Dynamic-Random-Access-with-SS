function error = fun_error(auset,rel,est_channels,est_data,uset,channels,data,JJ)
error = struct();
% This function calculates the user detection error, root mean square error
% of channel estimates, symbol error rate
rel_UE = auset(rel);
%% AER Estimation
fa = sum(~ismember(auset,uset)); %false alarm
md = sum(~ismember(uset,auset)); %misdetection
aer = (fa + md)/length(uset);
error.aer = aer;

%% NNMSE (Net Normalized Mean Square Error) Estimation
nnmse = norm(est_channels.'-channels(rel_UE))^2/norm(channels(rel_UE))^2;
error.nnmse = nnmse;
%% SER Estimation
[c,ia,ib] = intersect(uset,rel_UE); % Both reliable UE and actually transmitted
xor1 = xor(real(data(ia,:)),real(est_data(ib,:)));
xor2 = xor(imag(data(ia,:)),imag(est_data(ib,:)));
final = xor(xor1,xor2);
total_error = sum(final(:)) + 0*JJ*fa + 0*JJ*md; % accounted error due to misdetection and false alarm
ser = total_error/(size(est_data,1)*JJ);
error.ser = ser;
end