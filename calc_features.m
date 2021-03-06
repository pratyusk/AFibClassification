%% uses getFeatures.m to calculate
%go through all files and run getFeatures to preprocess and extract
%features
%---> [PHI,Y] matrix with features,
% calculate featurs and save them
clear all;
close all;
clc;
load('reference_test_lab_A_N.mat');
data = reference_test_lab;
Phi = [];


%% constants for feature calculation
plottf = 0; % no plots!
Fs = 300;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing design
    %% BP
    f_cutoff_low = 0.5; %to remove base 
    
    f_cutoff_high = 40; % in case there is
    
    
    order = 2;
    [b1,a1] = butter(order,[f_cutoff_low,f_cutoff_high]*2/Fs,'bandpass'); % Butterworth filter (low-pass by defualt). The first input is the filter oreder, the 2nd one is the cutoff frequency in the form (freq in Hz)*2/Fs
    %% double meadian filter
    %200msec = 60 samples; 600msec 180 samples
    n_median_200 = 60;
    n_median_600 = 180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% calculation of features
    
for i= 1:length(data)
   %calc featurs for each entry 
   filename = char(data(i,1));
   [tm,signal,Fs,siginfo]=rdmat(filename); % read from .mat file
   
   %%%%%%%%%%%%%%%%%PREPROCESS
%%%%%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%%%%%
%% remove mean
signal = signal - mean(signal);

%% apply butterworth BP filter designed in calc_features.m
signal = filter(b1,a1,signal); %apply BP filter

%% double meadian filter
%200msec = 60 samples; 600msec 180 sample
double_median = medfilt1(medfilt1(signal,n_median_200),n_median_600); %1D double median filter 200msec then 600msec
signal = signal - double_median;
   
 

    %calculate P,QRS features
    [R_idx,R_on_idx,P_idx,P_on_idx,P_off_idx,ann,type,num,RR_int_mean, RR_int_std, ratio_P2R,P_R_dist_mean,P_R_dist_std,P_R_on_dist_mean,P_R_on_dist_std] = ... 
    getFeatures_RR(filename,tm,signal,Fs,siginfo, plottf); %RR features
    
  
   %save data from ecgpuwave/rdann/gqrs
   ecgpuwave_data(i).ann = ann; ecgpuwave_data(i).num = num; ecgpuwave_data(i).type = type;
   ecgpuwave_data(i).class = cell2mat(data(i,2)); ecgpuwave_data(i).label = filename;
   ecgpuwave_data(i).signal = signal;
   
   
   Phi(i,1) = cell2mat(data(i,2)); %label
   
   
   % 9 PQRS features
        Phi(i,2:8) = ...
       [RR_int_mean, RR_int_std, ratio_P2R,P_R_dist_mean,P_R_dist_std,P_R_on_dist_mean,P_R_on_dist_std]; %RR features
    
  
   
   i
                                                %P-wave
   
       
end
                                                                                                                                                                                                         