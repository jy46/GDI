% Initialize
clear all
close all

% Add relevant paths
addpath ccdi_mat
addpath CTM_DI_package/Supporting_functions/
addpath CTM_DI_package/CTM_DI/

% User parameters
time_step   = 0.001;
M           = 3; % memory parameter
bin_width   = 35;  
boot_iter   = 2; % number of bootstrap iterations to use for classifier
time_max    = 300; % Maximum time to extract in seconds

% Load data
SNNAP_output = load('arb.out');
SNNAP_output(:,end) = []; % Output had a duplicate neuron, so remove
SNNAP_output(find(round(SNNAP_output(:,1),4)==time_max)+1:end,:) = []; % Keep to time_max

% Extract data
time = SNNAP_output(:,1);
data = SNNAP_output(:,2:end);

% Plot
for ii=1:size(data,2)
    plot(time, ii-0.5-0.5*(data(:,ii)/max(abs(data(:,ii)))))
    hold on
end
hold off
set(gca,'YDir','reverse')
yticks([0.5:1:(size(data,2)-0.5)])
ylabel('Neuron')
xlabel('Time (s)')
ylim([0 size(data,2)+0.25])
title('Raw Waveforms')

% Extract spike times
for ii=1:size(data,2)
    [peaks{ii},spike_times{ii}] = ...
        findpeaks(data(:,ii), time,'MinPeakHeight',0); % >0
    num_peaks(ii) = length(spike_times{ii});
end

% Plot spikes overlaying voltage waveforms
figure
for ii=1:size(data,2)
    plot(time, ii-0.5-0.5*(data(:,ii)/max(abs(data(:,ii)))))
    hold on
    scatter(spike_times{ii},ii-0.5-0.5*(peaks{ii}/max(abs(data(:,ii)))),'k*')
end
hold off
set(gca,'YDir','reverse')
yticks([0.5:1:(size(data,2)-0.5)])
ylabel('Neuron')
xlabel('Time (s)')
ylim([0 size(data,2)+0.25])
title('Spike Detection')

% TRUE RASTER PLOT !
figure
for ii=1:size(data,2)
    line([spike_times{ii} spike_times{ii}]', ...
         ii-0.5+[zeros(length(spike_times{ii}),1)+0.1 ones(length(spike_times{ii}),1)-0.1]',...
         'Color','k')
    hold on
end

hold off
yticks(1:11)
ylabel('Neuron')
xlabel('Time (s)')
ylim([0.5 11.5])
xlim([50 60])
title('Raster of Spike Times')


% Bin spikes
spike_times_binned = zeros(size(data));
for ii=1:length(spike_times)
    spike_times_binned(round(spike_times{ii}/time_step),ii) = 1;
end

% Plot binned spikes, also plot spike times from earlier for verification
% that binning worked
figure
for ii=1:size(data,2)
    plot(time, ii-0.5-0.9*spike_times_binned(:,ii))
    hold on
    scatter(spike_times{ii},ii-0.5-0.9*ones(size(spike_times{ii})),'k*')
end
hold off
set(gca,'YDir','reverse')
yticks([0:1:(size(data,2)-1)])
ylabel('Neuron')
xlabel('Time (s)')
ylim([min(ylim)-0.25 max(ylim)+0.25])
title('Raster of Binned Spikes')

% Further bin spike times
for ii=1:size(spike_times_binned,2)
    spike_times_binned_further(:,ii) = ...
        histcounts(find(spike_times_binned(:,ii)),...
            1:bin_width:length(spike_times_binned(:,ii)))';
end

% Plot further binned spike times
figure
time_further_binned = time(1:bin_width:end);
time_further_binned(end) = []; % Last bin edge not used
for ii=1:size(spike_times_binned_further,2)
    plot(time_further_binned,ii+spike_times_binned_further(:,ii))
    hold on
end
hold off
title('Raster of Further Binned Spikes')

% Estimate CTM-DI
% First convert spike_times 1) from a cell to zero padded
% array and 2) from sec to ms
spike_times_ms_array = zeros(length(time),length(spike_times));
for ii=1:length(spike_times)
    spike_times_ms_array(1:length(spike_times{ii}),ii) = ...
        spike_times{ii}*(1e3);
end
[CMatrix, w_i, HMatrix] = CTM_spiketime_wrapper(spike_times_ms_array,M,bin_width,1);

% Sign inference
[connection_sign, connection_sign_regular] = sign_inference(spike_times_binned_further,M);

% Estimate CCDI
[DI, DI_list] = di_compute(spike_times_binned_further,M,0,boot_iter);

DI(DI<0) = 0;

DI_norm = DI*(1/log(2))./HMatrix;
DI_norm(DI_norm<0.01) = 0;
DI_thresholded = (DI_norm/(1/log(2))).*HMatrix;

thresh = eps;

DI_cond_post = di_compute_post(DI_thresholded,thresh,M,spike_times_binned_further,boot_iter);
DI_cond_post(DI_cond_post<0)=0;
DI_cond_post_norm = DI_cond_post*(1/log(2))./HMatrix;

DI_cond_post_norm(logical(eye(11))) = 0;
DI_norm(logical(eye(11))) = 0;

%% TRUE CONNECTIVITY
true_connectivity = zeros(size(spike_times_binned_further,2),...
                          size(spike_times_binned_further,2));
true_connectivity(6,[1 3 5]) = [1 1 -1];
true_connectivity([6 9],2) = 1;
true_connectivity(8,2) = -1;
true_connectivity([3 5],10) = 1;
true_connectivity(1,4) = 1;
true_connectivity(4,11) = 1;

%% PLOT RESULTS
figure

subplot(1,3,1)
imagesc(true_connectivity)
colormap cool
colorbar
title('True Connectivity')

subplot(1,3,2)
imagesc(DI_norm)
colormap cool
colorbar
c_save(2,:) = caxis;
title('Estimated DI')

subplot(1,3,3)
imagesc(DI_cond_post_norm)
colormap cool
colorbar
c_save(3,:) = caxis;
title('Estimated GDI')

for ii=2:3
    subplot(1,3,ii)
    caxis([-max(abs(c_save(:))) max(abs(c_save(:)))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
