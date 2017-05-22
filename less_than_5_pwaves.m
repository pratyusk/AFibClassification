load('ecgpuwave_train_data_debug.mat')
load('Phi_debug.mat')
for i=1:length(Phi)
    
    P_wave_num(i) = length(find(~isnan(ecgpuwave_data(i).P_idx)));
    P_on_wave_num(i) = length(find(~isnan(ecgpuwave_data(i).P_on_idx)));
    P_off_wave_num(i) = length(find(~isnan(ecgpuwave_data(i).P_off_idx)));

end
% labels = Phi(:,1);
% boxplot( P_wave_num,labels)
length(find(Phi(find(P_wave_num < 5),1) ==0))
length(find(Phi(find(P_wave_num < 5),1) ==1))
length(find(Phi(find(P_on_wave_num < 5),1) ==0))
length(find(Phi(find(P_on_wave_num < 5),1) ==1))
length(find(Phi(find(P_off_wave_num < 5),1) ==0))
length(find(Phi(find(P_off_wave_num < 5),1) ==1))