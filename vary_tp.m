clc;
close all;
clear;
% Input order
% N,M,J,p,tp,radius,st,JJ
rng('shuffle')
%% Defining simulation parameter
N = 128;
M = 64;
Mlamb = [-28 -28 -28 -27 -27 -27 -26]; % thresholds for varying tp
J = 9;
p = 0.1;
tp = 3:3:21;
radius = 0.2;
st = 0.01;
JJ = 5;
pfa = 1e-5;
constellation = [1 1j -1 -1j]; %QPSK

AER_ESP = [];
NNMSE_ESP = [];
SER_ESP = [];
AER_SPICE = [];
NNMSE_SPICE = [];
SER_SPICE = [];
    rel_esprit = 0;
    rel_spice = 0;
for num=1:200000
    aer_esp = zeros(length(tp),1);
    nnmse_esp = zeros(length(tp),1);
    nser_esp = zeros(length(tp),1);
    aer_spice = zeros(length(tp),1);
    nnmse_spice = zeros(length(tp),1);
    nser_spice = zeros(length(tp),1);

    for mm = 1:length(tp)
        %% Generating spreading sequence
        phi = exp(1i*2*pi*(0:M-1)'/N*(0:N-1));
        phi = phi*diag(1./vecnorm(phi));
        %% Generating data
        params = mtc_data(N,M,J,p,tp(mm),radius,st,JJ,constellation);
        true_channels = params.channels;
        true_data = params.data;
        %% ESPRIT based Dynamic Random Access
        tstart = tic;
        [auset_esp,varest] = esprit_aud(params.snaps,N);
        [auset_ref_esp,act_mat_esp,XX_esp] = act_detect(auset_esp,phi,params.yn,pfa,varest);
        if sum(sum(act_mat_esp,2)==0)>0
            1;
        end
        [est_channels_esp,rel_esp] = channel_estimator(XX_esp,act_mat_esp,1);
        data_esp = data_detection(XX_esp(rel_esp,:),act_mat_esp(rel_esp,:),est_channels_esp,constellation); %detecting data of reliable UEs
        %data_esp = data;
        esprit_error = fun_error(auset_ref_esp,rel_esp,est_channels_esp,data_esp,params.uset,true_channels,true_data,J);
        rel_esprit = rel_esprit + length(rel_esp);
        %% Updating error vector for ESPRIT
        aer_esp(mm) = esprit_error.aer;
        nnmse_esp(mm) = esprit_error.nnmse;
        nser_esp(mm) = esprit_error.ser;
%        auset_esp = auset_ref;
        
        
        
        %% SPICE based Dynamic Random Access
        tstart = tic;
        auset_spice = spice_aud(params.snaps,phi,Mlamb(mm));
        [auset_ref,act_mat,XX] = act_detect(auset_spice,phi,params.yn,pfa,varest);
        if sum(sum(act_mat,2)==0)>0
            1;
        end
        [est_channels,rel] = channel_estimator(XX,act_mat,1);
        data = data_detection(XX(rel,:),act_mat(rel,:),est_channels,constellation); %detecting data of reliable UEs
        data_spice = data;
        spice_error = fun_error(auset_ref,rel,est_channels,data,params.uset,true_channels,true_data,J);
        rel_spice = rel_spice + length(rel);
        
        %% Updating error vector for SPICE
        aer_spice(mm) = spice_error.aer;
        nnmse_spice(mm) = spice_error.nnmse;
        nser_spice(mm) = spice_error.ser;
    end
    AER_ESP = [AER_ESP aer_esp];
    NNMSE_ESP = [NNMSE_ESP nnmse_esp];
    SER_ESP = [SER_ESP nser_esp];
    AER_SPICE = [AER_ESP aer_spice];
    NNMSE_SPICE = [NNMSE_SPICE nnmse_spice];
    SER_SPICE = [SER_SPICE nser_spice];
end
save('vary_tp.mat')
%% Producing meaningful statistics by weeding out NaN and outliers

mean(AER_ESP,2)
mean(AER_SPICE,2)
mean(stat_refine(NNMSE_ESP,0.05),2)
mean(stat_refine(NNMSE_SPICE,0.05),2)
mean(stat_refine(SER_ESP,0.05),2)
mean(stat_refine(SER_SPICE,0.05),2)
rel_esprit
rel_spice