%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPG TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
clear all
close all

% Add relevant paths
addpath ccdi_mat CTM_DI_package/CTM_DI CTM_DI_package/Supporting_functions

% User parameters
M           = 3;
bin_width   = 50; % Bin width in ms
boot_iter   = 5;
time_max    = 400; % Maximum time to extract in seconds
network_sel = 3; 
file_path   = 'spiketimes2.mat';

% Load data
SNNAP_output = load(file_path);

% Extract data
time_step      = SNNAP_output.srate; % sampling rate in seconds

rows_of_spikes_to_remove = find((time_step*SNNAP_output.spikes(:,1))>time_max); % find >time_max
SNNAP_output.spikes(rows_of_spikes_to_remove,:) = []; % remove >time_max

all_sim_spikes = SNNAP_output.spikes; % raw data from all simulations
all_sim_spikes(:,2) = all_sim_spikes(:,2)-min(all_sim_spikes(:,2))+1; % Ensure first neuron is index 1
selected_sim_spikes_indices = find(all_sim_spikes(:,3)==network_sel); % indices for sel sim
selected_sim_spikes = all_sim_spikes(selected_sim_spikes_indices,1:2);

% Convert selected_sim_spikes from time indices to time in seconds
num_neurons = length(unique(selected_sim_spikes(:,2)));
time = time_step*(1:max(selected_sim_spikes(:,1))); % Time in seconds
spike_times_binned = zeros([length(time) num_neurons]);
for ii=1:num_neurons
    current_neuron_spike_indices = find(selected_sim_spikes(:,2)==ii);
    current_neuron_spike_time_indices = selected_sim_spikes(current_neuron_spike_indices,1);
    spike_times_binned(current_neuron_spike_time_indices,ii) = 1;
    spike_times_ms_array(1:length(current_neuron_spike_time_indices),ii) = ...
        (1e3)*time(current_neuron_spike_time_indices); % conv to ms
end % time (which is in sec) corresponds to spike_times_binned, i.e. not really binning

% Plot spike times 
figure
for ii=1:num_neurons
    current_spike_times = find(spike_times_binned(:,ii));
    current_spike_times(time(current_spike_times)>time_max)=[];
    line([time(current_spike_times)' time(current_spike_times)']',...
       ii-1+[0.1+zeros(length(current_spike_times),1) ones(length(current_spike_times),1)-0.1]',...
       'Color','k')
    hold on
end
hold off
yticks((1:num_neurons)-0.5)
yticklabels(SNNAP_output.nname)
ylim([0 num_neurons])
xlabel('Time (s)')
ylabel('Neuron')
title('CPG Raw Raster Plot')

% Further bin spike times
for ii=1:num_neurons
    spike_times_binned_further(:,ii) = ...
        histcounts(find(spike_times_binned(:,ii)),...
            1:round((bin_width*0.001)/time_step):length(time))';
end
time_binned_further = time(1:round((bin_width*0.001)/time_step):end);

% Plot further binned spike times
figure
for ii=1:num_neurons
    current_spike_times = find(spike_times_binned_further(:,ii));
    current_spike_times(time_binned_further(current_spike_times)>time_max)=[];
    line([time_binned_further(current_spike_times)' time_binned_further(current_spike_times)']',...
       ii-1+[0.1+zeros(length(current_spike_times),1) ones(length(current_spike_times),1)-0.1]',...
       'Color','k')
    hold on
end
hold off
yticks((1:num_neurons)-0.5)
yticklabels(SNNAP_output.nname)
ylim([0 num_neurons])
xlabel('Time (s)')
ylabel('Neuron')
title('CPG Binned Further Raster Plot')

% Estimate CTM-DI
[CMatrix, w_i, HMatrix] = CTM_spiketime_wrapper(spike_times_ms_array,M,bin_width,1);

% Sign inference
[connection_sign, connection_sign_regular] = ...
    sign_inference(spike_times_binned_further,M);

% Estimate CCDI
warning('WARNING: Setting entries with >1 spike to 1')
spike_times_binned_further(spike_times_binned_further>1) = 1;

[DI, DI_list] = di_compute(spike_times_binned_further,M,0,boot_iter);

DI(DI<0) = 0;

DI_norm = DI*(1/log(2))./HMatrix;
DI_norm(DI_norm<0.01) = 0;
DI_thresholded = (DI_norm/(1/log(2))).*HMatrix;

thresh = eps;

DI_cond_post = di_compute_post(DI_thresholded,thresh,M,spike_times_binned_further,boot_iter);
DI_cond_post(DI_cond_post<0)=0;
DI_cond_post_norm = DI_cond_post*(1/log(2))./HMatrix;

DI_cond_post_norm(logical(eye(13))) = 0;
DI_norm(logical(eye(13))) = 0;

%% LOAD TRUE CONNECTIVITY
load('Actual_90820_2.mat')

%% PLOT RESULTS
figure

subplot(1,3,1)
imagesc(C)
colormap cool
colorbar
caxis([-max(abs(caxis)) max(abs(caxis))])
title('True Connectivity')

subplot(1,3,2)
imagesc(connection_sign_regular.*DI_norm)
colormap cool
colorbar
c_save(2,:) = caxis;
title('Estimated DI')

subplot(1,3,3)
imagesc(connection_sign.*DI_cond_post_norm)
colormap cool
colorbar
c_save(3,:) = caxis;
title('Estimated GDI')

for ii=2:3
    subplot(1,3,ii)
    caxis([-max(abs(c_save(:))) max(abs(c_save(:)))])
end

for ii=1:3
    subplot(1,3,ii)
    xticks(1:13)
    yticks(1:13)
    xticklabels(nname)
    yticklabels(nname)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
