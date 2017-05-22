%% uses getFeatures_p_sim.m to calculate
%---> [PHI,Y] matrix with features,
% calculate featurs and save them
% The features are: 
% mean and standard deviation of similarity of morphology of p waves (mean probably better)

clear all;
close all;
clc;
% load('ecgpuwave_train_data_debug.mat');
load('test_after_debug.mat'); %Phi matrix
start_idx = 1;

%% constants for feature calculation
plottf = 0; % no plots!
scale_sim = 10000; % for p wave similarity feature enhacement
sim_mean_base = 0.04; % as a base for cases if no p waves based on a few 
Fs = 300;
%trials of AFib examples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% calculation of features
    
for i=start_idx:length(ecgpuwave_data)
% for i=436:436
   %calc featurs for each entry 
   filename = ecgpuwave_data(i).label;
   signal = ecgpuwave_data(i).signal;
   tm = [0:1/Fs:(length(signal)-1)*(1/Fs)];
   R_idx = ecgpuwave_data(i).R_idx;
   R_on_idx = ecgpuwave_data(i).R_on_idx;
   P_idx = ecgpuwave_data(i).P_idx;
   P_on_idx = ecgpuwave_data(i).P_on_idx;
   P_off_idx = ecgpuwave_data(i).P_off_idx;
   
 
    %calculate P,QRS features
   [sim_P_wave_mean,sim_P_wave_std,sim_P_wave_mean_tri,sim_P_wave_std_tri] = ... 
       getFeatures_p_sim(filename,tm,signal,scale_sim,sim_mean_base,R_idx,R_on_idx,P_idx,P_on_idx,P_off_idx, plottf); %P sim calc
%    getFeatures_p_sim(filename,tm,signal,scale_sim,sim_mean_base,R_idx,R_on_idx,P_idx,P_on_idx,P_off_idx, plottf);
   Phi(i,9:10) = [sim_P_wave_mean,sim_P_wave_std];
   Phi(i,11:12) = [sim_P_wave_mean_tri,sim_P_wave_std_tri];
%    %calculate frequency bins features
%    [power_bins,energy_bins] = getFeatures_f_bins(filename,tm,signal,Fs,siginfo,binRange); %freq bins features
%    
% %    save data from ecgpuwave/rdann/gqrs
%    ecgpuwave_data(i).ann = ann; ecgpuwave_data(i).num = num; ecgpuwave_data(i).type = type; 
%    
%    
%    Phi(i,1) = cell2mat(data(i,2)); %label
%    
% %    frequ bin features
%    Phi(i,2:2+length(power_bins)+length(energy_bins)-1) = [power_bins,energy_bins];
%    
% %    9 PQRS features
%    Phi(i,2+length(power_bins)+length(energy_bins):2+length(power_bins)+length(energy_bins)+8) = ...
%        [RR_int_mean, RR_int_std, ratio_P2R,P_R_dist_mean,P_R_dist_std,P_R_on_dist_mean,P_R_on_dist_std,sim_P_wave_mean,sim_P_wave_std]; %RR features
%     
%    i
                                                %P-wave
   
       
end
    

save('Phi_test_with_sim_P.mat','Phi');