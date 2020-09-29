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

% Make plots tabs in figure window
set(0,'DefaultFigureWindowStyle','docked')

% User parameters
M           = 3;
bin_width   = 75; % Bin width in ms
boot_iter   = 100;
time_max    = 400; % Maximum time to extract in seconds
%thresh      = 0.005; % Threshold to apply to DI_thresholded and for condit
network_sel = 3; % sim 1: 20%, 2 is 15% and 3 is 10% noise
file_path   = 'CPG/spiketimes2.mat';

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

% For sake of old code, make spike_times_binned_further=smoothed_spike_times
smoothed_spike_times = spike_times_binned_further;

% Estimate CTM-DI
[CMatrix, w_i, HMatrix] = CTM_spiketime_wrapper(spike_times_ms_array,M,bin_width,1);

% Compute + plot partial cross correlation
r_partial = nan(size(smoothed_spike_times,2),size(smoothed_spike_times,2),(2*M)+1);
connection_sign = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_sum = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_regular = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_regular_sum = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
figure
load('Actual_90820_2.mat')
for ii=1:size(smoothed_spike_times,2)
    for jj=1:size(smoothed_spike_times,2)
        if ii~=jj
            [r_partial, r_regular] = partial_xcorr(smoothed_spike_times, ii, jj, M);
            r_partial_causal = r_partial(1:M);
            r_regular_causal = r_regular(1:M);
            
            if C(ii,jj)~=0
                subplot(size(smoothed_spike_times,2),...
                        size(smoothed_spike_times,2),...
                        jj+((ii-1)*size(smoothed_spike_times,2)))
                plot(-M:-1,r_partial_causal)
                xlim([-M -1])
            end
            
            I = find(max(abs(r_partial_causal))==abs(r_partial_causal));
            connection_sign(ii,jj) = sign(r_partial_causal(I));
            connection_sign_sum(ii,jj) = sum(r_partial_causal);
            
            I = find(max(abs(r_regular_causal))==abs(r_regular_causal));
            if length(I)>1
                warning(sprintf('Joe: %i %i Regular corr has multiple peaks',ii,jj));
                I = I(1);
            end
            connection_sign_regular(ii,jj) = sign(r_regular_causal(I));
            connection_sign_regular_sum(ii,jj) = sum(r_regular_causal);            
        end
    end
end

% Estimate CCDI
warning('WARNING: Setting entries with >1 spike to 1')
smoothed_spike_times(smoothed_spike_times>1) = 1;

[DI, DI_list] = di_compute(smoothed_spike_times,M,0,boot_iter);

DI(DI<0) = 0;

DI_norm = DI*(1/log(2))./HMatrix;
DI_norm(DI_norm<0.01) = 0;
DI_thresholded = (DI_norm/(1/log(2))).*HMatrix;

thresh = eps;

DI_cond_post = di_compute_post(DI_thresholded,thresh,M,smoothed_spike_times,boot_iter);
DI_cond_post(DI_cond_post<0)=0;
DI_cond_post_norm = DI_cond_post*(1/log(2))./HMatrix;

DI_cond_post_norm(logical(eye(13))) = 0;
DI_norm(logical(eye(13))) = 0;

% SORTING
[DI_norm_SORTED, I_DI_norm_SORTED] = ...
    sort(DI_norm(:));
[DI_cond_post_norm_SORTED, I_DI_cond_post_norm_SORTED] = ...
    sort(DI_cond_post_norm(:));

%% GRAPH

% FIll data structure graph_data to iterate over for plotting graphs
% graph_data{1} = C;
% graph_data{2} = DI_norm;
% graph_data{3} = DI_cond_post_norm;
% 
% % Loop through graph_data and plot
% figure
% for ii=1:3
%     [s, t] = find(graph_data{ii}~=0);
%     clearvars weights
%     for jj=1:length(s)
%         weights(jj) = abs(graph_data{ii}(s(jj),t(jj)));
%     end
%     
%     G = digraph(s,t,weights);
%     LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
%     
%     subplot(1,3,ii)
%     if ii==1
%         P(ii) = plot(G,'LineWidth',LWidths);
%     else
%         P(ii) = plot(G,'LineWidth',LWidths,'XData',P(1).XData,'YData',P(1).YData);
%     end
%     
%     P(ii).MarkerSize = 8;
%     
% end


% DI(logical(eye(length(DI)))) = 0;
% DI = connection_sign.*DI;
% DI

save(sprintf('CPG/runs_saved/M%i_BW%i_%is_Net%i_norm_B%i',...
             M,bin_width,time_max,network_sel,boot_iter))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%