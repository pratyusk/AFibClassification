function [ann,type,num,RR_int_mean, RR_int_std, ratio_P2R,P_R_dist_mean,P_R_dist_std,P_R_on_dist_mean,P_R_on_dist_std,sim_P_wave_mean,sim_P_wave_std] = ...
    getFeatures_RR_debug(ann, type, num,scale_sim,sim_mean_base, plottf); 

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



if nargin < 6
    plottf = 0;
end





%% %%%%%%%%%%%%%%%%%%%%%%%%clac QRS based features%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filenamedat = strcat(filename, '.dat'); 
% %if exist(filenamedat, 'file') == 0 %if not already existing
%     mat2wfdb(signal,strcat(filename,'wfdb'),Fs,[],siginfo.Units); %%change filename read .mat file, create 
%     %filename.dat and filenamewfdb.hea file for wfdb toolbox usage
    
    %% gqrs
%      gqrs(filename);
    %gqrs(strcat(filename,'wfdb')); %% change filename create qrs complex annotation
%     pause(5)
    %% ecgpuwave
%     ecgpuwave(filename,'ann',[],[],'qrs'); % create whole signal annotation

    %ecgpuwave(strcat(filename,'wfdb'),'ann',[],[],'qrs'); % create whole signal annotation
    %% cannot use both gqrs and sqrs at the same time since both create .qrs annotation
%end
%nwaves_p = rdann_new(strcat(filename,'wfdb'),'qrs'); %% change filename read annotation, output according to input type
% [ann,type,subtype,chan,num,comments] = rdann_new(filename,'ann',[],[],[]); %last input 'N'
%[ann,type,subtype,chan,num,comments] = rdann_new(strcat(filename,'wfdb'),'ann',[],[],[]); %last input 'N'
 % ann contains the real idx for each annotation point
 % type defines what kind of annotation contains symbols "(" (onset belonging to see num vector), 
% ")" (offset), "p" (p-peak), "N"(R peak)
%comment, chan...
% num tells to which wavetype an onset/offset corresponds, only makes sense
% if in type "(" or ")"


%% extract relevant points
% type vector contains symbols "(" (onset belonging to see num vector), 
% ")" (offset), "p" (p-peak), "N"(R peak) --> look in ann for the "real" idx
% of the event in type vectore 

%% R_idx

% R peaks
%idx in type annotation vector
R_in_type_idx = find(type == 'N');
%correspondind idx in signal
R_idx = ann(R_in_type_idx);



% QRS onset
R_on_idx = NaN(length(R_idx),1); % same length as R peak idx vector
%NaN if no onset for the R peak

% first element
if R_in_type_idx(1)-1 >0 % in case R is the first annotation
    
    if type(R_in_type_idx(1)-1) == '(' && num(R_in_type_idx(1)-1) == 1
        %if there is an onset annotation and it is labeled QRS for the 
        %first R wave
       R_on_idx(1) = ann(R_in_type_idx(1)-1); 
    end
    
end
%further elements
if length(R_idx) >1 % if we have more R waves than 1
for i = 2:length(R_idx)   
        if type(R_in_type_idx(i)-1) == '(' && num(R_in_type_idx(i)-1) == 1; 
            %is there an onset to the R peak?
            R_on_idx(i) = ann(R_in_type_idx(i)-1); 
        end
end
end


%% p_idx

%vectors same length as R_peak_idx
% peaks
P_idx = NaN(length(R_idx),1);
% onset p
P_on_idx = NaN(length(R_idx),1);
% offset p
P_off_idx = NaN(length(R_idx),1);

P_in_type_first = find(type(1:R_in_type_idx(1)) == 'p'); % firsr P wave before first R wave
   if ~isempty(P_in_type_first) %if there was a p wave before the first R wave
        P_idx(1) = ann(min(P_in_type_first));  %just in rare case  that 2 p in a row and max/min abitary
        idx_P_in_first = min(P_in_type_first);
        if idx_P_in_first-1 > 0 % case that p wave is the first sample
            if type(idx_P_in_first-1) == '(' && num(idx_P_in_first-1) == 0; 
                %is there an onset to the p wave and is it labeled P?
                 P_on_idx(1) = ann(idx_P_in_first-1); %then save the corresponding signal idx
            end
        end
        
        if type(idx_P_in_first+1) == ')' && num(idx_P_in_first+1) == 0; 
            %is there an offset to the P wave and is it labled P?
            P_off_idx(1) = ann(idx_P_in_first+1);  
        end
   
   
   
    end

% we look between 2 R peaks if there was a p wave. assign NaN if not or
% assign the right index(minimal in case there are more) if there was a
% pwave
if length(R_idx) > 1 %if there are more than 1 R peak
for i = 1:length(R_idx)-1
    % search in type vector between 2 R peak events for P peak events
    P_in_type = find(type(R_in_type_idx(i):R_in_type_idx(i+1)) == 'p'); %type_vec between 2 R peaks i and i+1 so before i+1
    % this P wave corresponds to R peak i+1
    if ~isempty(P_in_type) %if there was a p wave between this 2 R peaks
        % offset because we look after R peak i idx and before R idx i+1
        P_idx(i+1) = ann(min(R_in_type_idx(i) + P_in_type)-1);  %just in rare case  that 2 p in a row and max/min abitary
        idx_P_in_type = min(R_in_type_idx(i) + P_in_type)-1;
        
        %is there an onset to the p wave? and assigned as P
        if type(idx_P_in_type-1) == '(' && num(idx_P_in_type-1) == 0;  
            P_on_idx(i+1) = ann(idx_P_in_type-1); 
        end
        
        if idx_P_in_type+1 <= length(type) %p wave last annotation made
            %is there an offset to the p wave? 
        if type(idx_P_in_type+1) == ')'&& num(idx_P_in_type+1) == 0;  
            P_off_idx(i+1) = ann(idx_P_in_type+1); 
        end
        end
        
    end
end
end

%% features based on the points calculated above

%% R-R interval features

% RR interval in terms of samples

if length(R_idx) >1 % if it is possible to calculate a distance
    RR_interval = zeros(length(R_idx) - 1, 1); % get the peak interval in seconds
    RR_interval = R_idx(2:end) - R_idx(1:end-1);


    RR_int_mean = mean(RR_interval);%mean of RR intervals in sec

    if length(R_idx) > 11 %only if we have more than 5 distances it makes sense
        %to compute std
        RR_int_std = std(RR_interval); %std of RR intervals in sec
    else
        RR_int_std = NaN;

    end
end

%% P-R interval and ratio


% ratio p waves detected to R peaks detected in % 
ratio_P2R = 100* (length(P_idx)-length(find(isnan(P_idx))))/length(R_idx);

% P-R distance based on peaks
% only consider if there is a P
% wave%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_R_dist = zeros(length(R_idx),1);
P_R_dist = R_idx(find(~isnan(P_idx))) - P_idx(find(~isnan(P_idx)));

P_R_dist_mean = mean(P_R_dist);

%if there are less then 5 distances/P waves then set std to NaN
if length(P_R_dist > 5)
    P_R_dist_std = std(P_R_dist);
else
    P_R_dist_std = NaN;
end

% P-R distance based on onsets of p and
% QRS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%only if possible so if to each QRS onset there is an P onset

P_R_on_dist = R_on_idx(find(~isnan(P_on_idx) & ~isnan(R_on_idx))) - P_on_idx(find(~isnan(P_on_idx) & ~isnan(R_on_idx)));

P_R_on_dist_mean = mean(P_R_on_dist);
% if less than 5 distances then set std to NaN
P_R_on_dist_std = std(P_R_on_dist);


% % similarity
% % we calculate the similarity between 2 concecutive p waves
% % only possible if consecutive ones exists
% 
% sim_P_waves = zeros(length(R_idx));
% 
% %interpolate p wave, p onset and  p offset for each beat
% %distance
% p_r_dist = R_idx - P_idx;
% p_on_r_dist = R_idx - P_on_idx;
% p_off_r_dist = R_idx - P_off_idx;
% 
% for i=1:length(idx_exists_p)-1
%     peak2on = 0;  %differnce for peak to onset distance
%     peak2off = 0; %differnce for peak to offset distance
%     on2off = 0;   %difference for onset to offset distance
%     line1 = 0;     % to check peak and onset exists for both consecutive pwaves
%     line2 = 0;     % to check peak and offset exist for both consecutive pwaves
%     if idx_exists_p(i+1) - idx_exists_p(i)== 1 %if consecutive p waves check similarity
%        if P_on_idx(idx_exists_p(i)) ~= 0 && P_on_idx(idx_exists_p(i+1)) ~= 0 %if on to peak possible
%         line1 = 1; % peak and offset exist for both consecutive pwaves  
%         a = sqrt((signal(P_on_idx(idx_exists_p(i)))- signal(P_idx(idx_exists_p(i))))^2 + (tm(P_on_idx(idx_exists_p(i)))-tm(P_idx(idx_exists_p(i))))^2);
%         b = sqrt((signal(P_on_idx(idx_exists_p(i+1)))- signal(P_idx(idx_exists_p(i+1))))^2 + (tm(P_on_idx(idx_exists_p(i+1)))-tm(P_idx(idx_exists_p(i+1))))^2);
%         peak2on = (a-b)^2;
%        end
%        
%        if P_off_idx(idx_exists_p(i)) ~= 0 && P_off_idx(idx_exists_p(i+1)) ~= 0 %if off to peak possible
%          line2 = 1;     %peak and offset exist for both consecutive pwaves   
%         a = sqrt((signal(P_off_idx(idx_exists_p(i)))- signal(P_idx(idx_exists_p(i))))^2 + (tm(P_off_idx(idx_exists_p(i)))-tm(P_idx(idx_exists_p(i))))^2);
%         b = sqrt((signal(P_off_idx(idx_exists_p(i+1)))- signal(P_idx(idx_exists_p(i+1))))^2 + (tm(P_off_idx(idx_exists_p(i+1)))-tm(P_idx(idx_exists_p(i+1))))^2);
%         peak2off = (a-b)^2;
%        end
%        
%        if line1 == 1 && line2 == 1 %if onset and offset in both
%             a = sqrt((signal(P_off_idx(idx_exists_p(i)))- signal(P_on_idx(idx_exists_p(i))))^2 + (tm(P_off_idx(idx_exists_p(i)))-tm(P_on_idx(idx_exists_p(i))))^2);
%             b = sqrt((signal(P_off_idx(idx_exists_p(i+1)))- signal(P_on_idx(idx_exists_p(i+1))))^2 + (tm(P_off_idx(idx_exists_p(i+1)))-tm(P_on_idx(idx_exists_p(i+1))))^2);
%            on2off = (a-b).^2;
%        end
%        dist_vec = [peak2on peak2off on2off];
%        %only consider the  lines exsisting in both consecutive p waves
%        % for comparison reason normalize difference by number of lines
%        % evaluated
%        if length(find(dist_vec)) >0
%             sim_P_waves(i) = sqrt(sum(dist_vec))/length(find(dist_vec));
%        end
%     end
% end
% 
% sim_P_wave_mean = mean(sim_P_waves)*scale_sim; %for enhacement
% sim_P_wave_std = std(sim_P_waves)*scale_sim;
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
