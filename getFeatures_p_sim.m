function [sim_P_wave_mean,sim_P_wave_std,sim_P_wave_mean_tri,sim_P_wave_std_tri] = ...
            getFeatures_p_sim(filename,tm,signal,scale_sim,sim_mean_base,R_idx,R_on_idx,P_idx,P_on_idx,P_off_idx, plottf); %RR features

% Usage: [beatstd, beatstdiqr, signalstd, signalstdiqr, histRatio, iqrtRange, beatmean] = getFeatures(filename, plottf, binRange)
% if plottf = 0, function will not plot (default 0)
% if bin range not specified, default binRange = [0:0.1:2]

% input: the ecg .mat filename as a string, sampling frequency, butterworth filter coeff, median_filter lenghts
%. The function will read the .mat 
% tm,signal,Fs,siginfo from rdmat
% file and output the .dat wfdb format file and .hea annotation file, then 
% generate estimate qrs annotation file (.qrs). 

% output: standard deviation of the peak interval, the
% interval mean, indices of peak location



if nargin < 10
    plottf = 0;
end



% %%%%%%%%%%%%%%%%%PREPROCESS
% %%%%%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%%%%%
% %% remove mean
% signal = signal - mean(signal);
% 
% %% apply butterworth BP filter designed in calc_features.m
% signal = filter(b1,a1,signal); %apply BP filter
% 
% %% double meadian filter
% %200msec = 60 samples; 600msec 180 sample
% double_median = medfilt1(medfilt1(signal,n_median_200),n_median_600); %1D double median filter 200msec then 600msec
% signal = signal - double_median;





%% %%%%%%%%%%%%%%%%%%%%%%%%clac QRS based features%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenamedat = strcat(filename, '.dat'); 
% %if exist(filenamedat, 'file') == 0 %if not already existing
%     mat2wfdb(signal,strcat(filename,'wfdb'),Fs,[],siginfo.Units); %%change filename read .mat file, create 
%     %filename.dat and filenamewfdb.hea file for wfdb toolbox usage
%     
%     %% gqrs
%      gqrs(filename);
%     %gqrs(strcat(filename,'wfdb')); %% change filename create qrs complex annotation
% %     pause(5)
%     %% ecgpuwave
%     ecgpuwave(filename,'ann',[],[],'qrs'); % create whole signal annotation
% 
%     %ecgpuwave(strcat(filename,'wfdb'),'ann',[],[],'qrs'); % create whole signal annotation
%     %% cannot use both gqrs and sqrs at the same time since both create .qrs annotation
% %end
% %nwaves_p = rdann_new(strcat(filename,'wfdb'),'qrs'); %% change filename read annotation, output according to input type
% [ann,type,subtype,chan,num,comments] = rdann_new(filename,'ann',[],[],[]); %last input 'N'

%[ann,type,subtype,chan,num,comments] = rdann_new(strcat(filename,'wfdb'),'ann',[],[],[]); %last input 'N'
 % ann contains the real idx for each annotation point
 % type defines what kind of annotation contains symbols "(" (onset belonging to see num vector), 
% ")" (offset), "p" (p-peak), "N"(R peak)
%comment, chan...
% num tells to which wavetype an onset/offset corresponds, only makes sense
% if in type "(" or ")"


%% similarity of P waves
% we calculate the similarity between 2 concecutive p waves
% only possible if consecutive ones exists

sim_P_waves = zeros(length(R_idx),1);
sim_P_waves(:) = sim_mean_base;

%interpolate p wave, p onset and  p offset for each beat
%distance
if (length(find(~isnan(P_on_idx))) >= 2) && (length(find(~isnan(P_off_idx))) >= 2)

    % interpolate p wave peaks, onset, offset
%     P_idx_interp = floor(interp1(find(~isnan(P_idx)), P_idx(~isnan(P_idx)), 1:length(P_idx)));
    P_on_idx_interp = floor(interp1(find(~isnan(P_on_idx)), P_on_idx(~isnan(P_on_idx)), 1:length(P_on_idx)));
    P_off_idx_interp = floor(interp1(find(~isnan(P_off_idx)), P_off_idx(~isnan(P_off_idx)), 1:length(P_off_idx)));


    % index of p where value is not NaN
    % idx_exists_p = find(~isnan(P_idx_interp));

    cross_correlation = zeros(length(P_on_idx_interp)-1,1);
    cross_correlation(:) = NaN;
    % (length(P_idx_interp)-1)
    for i=1:(length(P_on_idx_interp)-1)
        line1 = 0;     % to check peak and onset exists for both consecutive pwaves
        line2 = 0;     % to check peak and offset exist for both consecutive pwaves
        if (P_off_idx_interp(i) - P_on_idx_interp(i) > 0)
            signal_before = signal(P_on_idx_interp(i):P_off_idx_interp(i));
            line1 = 1;
        end
        if (P_off_idx_interp(i+1) - P_on_idx_interp(i+1) > 0)
            signal_after = signal(P_on_idx_interp(i+1):P_off_idx_interp(i+1));
            line2 = 1;
        end
        if (line1 == 1 && line2 == 1)
            cross_correlation(i) = max(abs(xcorr(signal_before, signal_after)));
        end
        i
        filename
    end

    sim_P_wave_mean = mean(cross_correlation(find(~isnan(cross_correlation)))); %*scale_sim; %for enhacement
    sim_P_wave_std = std(cross_correlation(find(~isnan(cross_correlation)))); %*scale_sim;
else 
    sim_P_wave_mean = NaN;
    sim_P_wave_std = NaN;
end



sim_P_waves_tri = zeros(length(R_idx),1);
sim_P_waves_tri(:) = sim_mean_base;
if (length(find(~isnan(P_idx))) >= 2) && (length(find(~isnan(P_on_idx))) >= 2) && (length(find(~isnan(P_off_idx))) >= 2)
    P_idx_interp = floor(interp1(find(~isnan(P_idx)), P_idx(~isnan(P_idx)), 1:length(P_idx)));
    P_on_idx_interp = floor(interp1(find(~isnan(P_on_idx)), P_on_idx(~isnan(P_on_idx)), 1:length(P_on_idx)));
    P_off_idx_interp = floor(interp1(find(~isnan(P_off_idx)), P_off_idx(~isnan(P_off_idx)), 1:length(P_off_idx)));


    idx_exists_p = find(~isnan(P_idx_interp));

    for i=1:length(idx_exists_p)-1
        peak2on = 0;  %differnce for peak to onset distance
        peak2off = 0; %differnce for peak to offset distance
        on2off = 0;   %difference for onset to offset distance
        line1 = 0;     % to check peak and onset exists for both consecutive pwaves
        line2 = 0;     % to check peak and offset exist for both consecutive pwaves
        if idx_exists_p(i+1) - idx_exists_p(i)== 1 %if consecutive p waves peaks, check similarity
           if ~isnan(P_on_idx_interp(idx_exists_p(i))) && ~isnan(P_on_idx_interp(idx_exists_p(i+1))) %if on to peak possible
            line1 = 1; % peak and offset exist for both consecutive pwaves  
            a = sqrt((signal(P_on_idx_interp(idx_exists_p(i)))- signal(P_idx_interp(idx_exists_p(i))))^2 + (tm(P_on_idx_interp(idx_exists_p(i)))-tm(P_idx_interp(idx_exists_p(i))))^2);
            b = sqrt((signal(P_on_idx_interp(idx_exists_p(i+1)))- signal(P_idx_interp(idx_exists_p(i+1))))^2 + (tm(P_on_idx_interp(idx_exists_p(i+1)))-tm(P_idx_interp(idx_exists_p(i+1))))^2);
            peak2on = (a-b)^2;
           end

           if ~isnan(P_off_idx_interp(idx_exists_p(i))) && ~isnan(P_off_idx_interp(idx_exists_p(i+1))) %if off to peak possible
             line2 = 1;     %peak and offset exist for both consecutive pwaves   
            a = sqrt((signal(P_off_idx_interp(idx_exists_p(i)))- signal(P_idx_interp(idx_exists_p(i))))^2 + (tm(P_off_idx_interp(idx_exists_p(i)))-tm(P_idx_interp(idx_exists_p(i))))^2);
            b = sqrt((signal(P_off_idx_interp(idx_exists_p(i+1)))- signal(P_idx_interp(idx_exists_p(i+1))))^2 + (tm(P_off_idx_interp(idx_exists_p(i+1)))-tm(P_idx_interp(idx_exists_p(i+1))))^2);
            peak2off = (a-b)^2;
           end

           if line1 == 1 && line2 == 1 %if onset and offset in both
                a = sqrt((signal(P_off_idx_interp(idx_exists_p(i)))- signal(P_on_idx_interp(idx_exists_p(i))))^2 + (tm(P_off_idx_interp(idx_exists_p(i)))-tm(P_on_idx_interp(idx_exists_p(i))))^2);
                b = sqrt((signal(P_off_idx_interp(idx_exists_p(i+1)))- signal(P_on_idx_interp(idx_exists_p(i+1))))^2 + (tm(P_off_idx_interp(idx_exists_p(i+1)))-tm(P_on_idx_interp(idx_exists_p(i+1))))^2);
               on2off = (a-b).^2;
           end
           dist_vec = [peak2on peak2off on2off];
           %only consider the  lines exsisting in both consecutive p waves
           % for comparison reason normalize difference by number of lines
           % evaluated
           if length(find(dist_vec)) >0
                sim_P_waves_tri(i) = sqrt(sum(dist_vec))/length(find(dist_vec));
           end
        end
    end

    sim_P_wave_mean_tri = mean(sim_P_waves_tri)*scale_sim; %for enhacement
    sim_P_wave_std_tri = std(sim_P_waves_tri)*scale_sim;
else
    sim_P_wave_mean_tri = NaN;
    sim_P_wave_std_tri = NaN;
end






% 
% 
% 
% 
% 
% %% plot signal and peaks detected
% if (plottf ~= 0)
%     figure;plot(tm,signal(:,1));hold on;grid on
% %     plot(tm(nwaves_p),signal(nwaves_p),'or')
% %     plot(tm(nwaves_R),signal(nwaves_R),'go')
%    plot(tm(R_idx),signal(R_idx),'go')
%     plot(tm(R_on_idx+1),signal(R_on_idx+1),'go')
%    
%     plot(tm(P_idx+1),signal(P_idx+1),'ro')
%     plot(tm(P_on_idx+1),signal(P_on_idx+1),'ro')
%     plot(tm(P_off_idx+1),signal(P_off_idx+1),'ro')
%     title(filename);
% end
% 
% 
% % %% calculate P wave features
% % nwave_p_interval = zeros(length(nwaves_p) - 1, 1); % get the peak interval in seconds
% % for i = 1:(length(nwaves_p) - 1)
% %     nwave_p_interval(i, 1) = tm(nwaves_p(i+1)) - tm(nwaves_p(i));
% % end
% % P_interval = nwave_p
% 
% 
% 
