function params = mtc_data(N,M,J,p,tp,radius,st,JJ,constellation)
params = struct();
AU = ceil(p*N);
ud = rand(N,1)*(radius-st) + st;
phi = exp(1i*2*pi*(0:M-1)'/N*(0:N-1));
phi = phi*diag(1./vecnorm(phi));
%phi = sqrt(10^(tp/10))*phi;
% Received signal generation
uset = sort(randsample(2:N,AU));
%uset = datasample(2:128,AU);
uphi = phi(:,uset); % matrix of codes used by UEs
rp = -128.1 - 36.7*log10(ud)+tp; % in dBm
rp = 10.^(rp/10); % in mW
pc = rp;
h = sqrt(1/2)*(normrnd(0,1,N,1) + 1i*normrnd(0,1,N,1)); % Rayleigh fading
h = h.*sqrt(rp); % Combining pathloss and channel fading
constellation = [1 1j -1 -1j]; %QPSK
% Designing activity matrix
A = zeros(AU,J);
AA = zeros(N,J);
data = zeros(AU,J);
ct = 1;
% Data is generated under the assumption that UE enters and exit within the
% random access opportunity 
for i = 1:length(uset)
    %find the starting of transmission
    io = randperm(J-JJ+1,1); % to make sure the UE enters and exit within RA opportunity
    %act_ind = mod(io:io+JJ-1,J) + 1; %active indices
    act_ind = io:io+JJ-1;
    %temp = constellation(randperm(length(constellation),JJ));
    temp = randsample(constellation,JJ,true);
    temp(1) = 1 ; %pilot
    %data = [data;temp];
    data(i,act_ind) = temp;
    A(ct,act_ind) = h(uset(i))*temp;
    AA(uset(i),act_ind) = h(uset(i))*temp;
    ct = ct + 1;
end

y = uphi*A;
yy = phi*AA; % Sparse representation problem
np = -110; %noise power in dBm for transmission bandwidth 1 MHz
noise = sqrt(1/2)*(normrnd(0,1,M,J) + 1i*normrnd(0,1,M,J))*(sqrt(10^(np/10))*eye(J,J)/2); % Not normalized by spreading sequence
yn = y + noise; % received signal in symbols
yyn = yy + noise;

%% Creating snapshots
switch M
    case 32
        m = 24;
    case 48
        m = 40;
    case 64
        m = 50;
    case 80
        m = 60;
    case 96
        m = 78;
    case 112
        m = 90;
    case 432
        m = 416;
end
snaps = [];
for i=1:J
    snaps = [snaps snapdim(yn(:,i),m)];
end
params.yn = yn;
params.ts = A; %data received in timeslot 
params.data = data; %data symbol sent by active users
params.snaps = snaps;
params.uset = uset;
params.channels = h;
end
