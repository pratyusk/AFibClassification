close all; clear all; clc;
load('ecgpuwave_train_data_debug.mat');
for k = 1:length(ecgpuwave_data)
type = ecgpuwave_data(k).type;
num = ecgpuwave_data(k).num;    
ann = ecgpuwave_data(k).ann;

% state P
P(k).state = NaN();

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
   P_on_in_type_first = find(type(1:R_in_type_idx(1)) == '('); %  onsets before first R wave
   P_off_in_type_first = find(type(1:R_in_type_idx(1)) == ')'); % offsets before first R wave
   
   % if there was on onset before the first R peak && there are onset(s) marked as P &&
   % and we didnt find one P onset yet
   if ~isempty(P_on_in_type_first) && length(find(num(P_on_in_type_first) == 0)) > 0 &&  isnan(P_on_idx(1))
       %then assign the latest P onset idx
       P_on_idx(1) = ann(max(find(num(P_on_in_type_first) == 0)));
   end
   
   % if there was on offset before the first R peak && there are offset(s) marked as P &&
   % and we didnt find one P offset yet
   if ~isempty(P_off_in_type_first) && length(find(num(P_off_in_type_first) == 0)) > 0 &&  isnan(P_off_idx(1))
       %then assign the latest P onset idx
       P_off_idx(1) = ann(max(find(num(P_off_in_type_first) == 0)));
   end
   
  
   
   %check state
   if ~isnan(P_idx(1)) && ~isnan(P_on_idx(1)) && ~isnan(P_off_idx(1)) 
       P(k).state(1) = 7;
        elseif ~isnan(P_idx(1)) && ~isnan(P_on_idx(1)) && isnan(P_off_idx(1))
        P(k).state(1) = 6;
       elseif ~isnan(P_idx(1)) && isnan(P_on_idx(1)) && ~isnan(P_off_idx(1))
         P(k).state(1) = 5;
        elseif ~isnan(P_idx(1)) && isnan(P_on_idx(1)) && isnan(P_off_idx(1))
         P(k).state(1) = 4;
        elseif isnan(P_idx(1)) && ~isnan(P_on_idx(1)) && ~isnan(P_off_idx(1))
         P(k).state(1) = 3;
        elseif isnan(P_idx(1)) && ~isnan(P_on_idx(1)) && isnan(P_off_idx(1))
         P(k).state(1) = 2;
        elseif isnan(P_idx(1)) && isnan(P_on_idx(1)) && ~isnan(P_off_idx(1))
         P(k).state(1) = 1;
   else 
         P(k).state(1) = 0;
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
    
   P_on_in_type = find(type(R_in_type_idx(i):R_in_type_idx(i+1)) == '('); %  onsets between 2 R peaks
   P_off_in_type = find(type(R_in_type_idx(i):R_in_type_idx(i+1)) == ')'); %  offsets between 2 R peaks
    
   
    % if there was one onsetbetween 2 R peaks && there are onset(s) marked as P &&
   % and we didnt find one P onset yet
   if ~isempty(P_on_in_type) && length(find(num(P_on_in_type) == 0)) > 0 &&  isnan(P_on_idx(i))
       %then assign the latest P onset idx
       P_on_idx(i) = ann(max(find(num(P_on_in_type) == 0)));
   end
   
   % if there was on offset between two R peaks && there are offset(s) marked as P &&
   % and we didnt find one P offset yet
   if ~isempty(P_off_in_type) && length(find(num(P_off_in_type) == 0)) > 0 &&  isnan(P_off_idx(i))
       %then assign the latest P onset idx
       P_off_idx(i) = ann(max(find(num(P_off_in_type) == 0)));
   end
   
  
   
   %check state
   if ~isnan(P_idx(i)) && ~isnan(P_on_idx(i)) && ~isnan(P_off_idx(i)) 
       P(k).state(i) = 7;
        elseif ~isnan(P_idx(i)) && ~isnan(P_on_idx(i)) && isnan(P_off_idx(i))
        P(k).state(i) = 6;
       elseif ~isnan(P_idx(i)) && isnan(P_on_idx(i)) && ~isnan(P_off_idx(i))
         P(k).state(i) = 5;
        elseif ~isnan(P_idx(i)) && isnan(P_on_idx(i)) && isnan(P_off_idx(i))
         P(k).state(i) = 4;
        elseif isnan(P_idx(i)) && ~isnan(P_on_idx(i)) && ~isnan(P_off_idx(i))
         P(k).state(i) = 3;
        elseif isnan(P_idx(i)) && ~isnan(P_on_idx(i)) && isnan(P_off_idx(i))
         P(k).state(i) = 2;
        elseif isnan(P_idx(i)) && isnan(P_on_idx(i)) && ~isnan(P_off_idx(i))
         P(k).state(i) = 1;
   else 
         P(k).state(i) = 0;
   end
    
end
end

end


% what happened ?
 state_P_0 = [];
 state_P_1 = [];

for i = 1:length(P)    
    if ecgpuwave_data(i).class == 'N'
   state_P_0  = [state_P_0, P(i).state]; %concanate   
    else   
   state_P_1  = [state_P_1, P(i).state]; %concanate 
    end
end

length(find(isnan(state_P_0)))
length(find(isnan(state_P_1)))

for i = 0:7
   i;
   disp(strcat('Normal: Number of samples with state : ',num2str(i),' = ',num2str(length(find(state_P_0 == i)))))
   disp(strcat('Normal: Percent of samples with state : ',num2str(i),' = ',num2str(length(find(state_P_0 == i))/length(state_P_0)*100)))
   disp(strcat('AFib: Number of samples with state : ',num2str(i),' = ',num2str(length(find(state_P_1 == i)))))
   disp(strcat('AFib: Percent of samples with state : ',num2str(i),' = ',num2str(length(find(state_P_1 == i))/length(state_P_1)*100)))
end