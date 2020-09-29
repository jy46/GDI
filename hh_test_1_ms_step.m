% Initialize
clear all
close all

% Add relevant paths
addpath hh_testing
addpath ccdi_mat
addpath CTM_DI_package/Supporting_functions/
addpath CTM_DI_package/CTM_DI/

% Make plots tabs in figure window
set(0,'DefaultFigureWindowStyle','docked')

% Colorbrewer
addpath DrosteEffect-BrewerMap-5b84f95/
set(0,'DefaultAxesColorOrder',brewermap(11,'Set1'))
colors = colororder;

% User parameters
time_step   = 0.001;
sigma       = 2;
M           = 3;
C           = 0;
bin_width   = 30;
boot_iter   = 100;
kernel_type = 'none';
time_max    = 200; % Maximum time to extract in seconds
network_sel = 'rockon';

% For Fig5 networks B, C
% bin_width   = 25;
% M           = 10;

% For Fig6 network _A1
% bin_width   = 60;
% M           = 5;

% For Fig6 network _A1, _A2, _A3
% bin_width   = 100;
% M           = 3;

% Load data
SNNAP_output = load(['hh_testing/network' network_sel '.smu.out']);
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
yticklabels({'A';'B';'C'})
ylabel('Neuron')
xlabel('Time (s)')
ylim([0 size(data,2)+0.25])

% Extract spike times
for ii=1:size(data,2)
    [peaks{ii},spike_times{ii}] = ...
        findpeaks(data(:,ii), time,'MinPeakHeight',0); % >0
    num_peaks(ii) = length(spike_times{ii});
end

% Plot spikes overlaying voltage waveforms
figure
for ii=1:size(data,2)
    plot(time, ii-0.5-0.5*(data(:,ii)/max(abs(data(:,ii)))), 'Color',...
         colors(ii,:))
    hold on
    scatter(spike_times{ii},ii-0.5-0.5*(peaks{ii}/max(abs(data(:,ii)))),'k*')
end
hold off
set(gca,'YDir','reverse')
yticks([0.5:1:(size(data,2)-0.5)])
yticklabels({'A';'B';'C'})
ylabel('Neuron')
xlabel('Time (s)')
ylim([0 size(data,2)+0.25])


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


% Bin spikes
spike_times_binned = zeros(size(data));
for ii=1:length(spike_times)
    spike_times_binned(round(spike_times{ii}/time_step),ii) = 1;
end

% Plot binned spikes, also plot spike times from earlier for verification
% that binning worked
figure
for ii=1:size(data,2)
    plot(time, ii-0.5-0.9*spike_times_binned(:,ii), 'Color',...
         colors(ii,:))
    hold on
    scatter(spike_times{ii},ii-0.5-0.9*ones(size(spike_times{ii})),'k*')
end
hold off
set(gca,'YDir','reverse')
yticks([0:1:(size(data,2)-1)])
yticklabels({'A';'B';'C'})
ylabel('Neuron')
xlabel('Time (s)')
ylim([min(ylim)-0.25 max(ylim)+0.25])

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

% Smooth binned spike times with Gaussian kernel
edges  = -3*sigma:1:3*sigma;
switch kernel_type
    case 'gaussian'
        kernel = normpdf(edges,0,sigma);
    case 'gamma'
        kernel = gampdf(edges,1,1);
    case 'none'
        kernel = zeros(1,length(edges));
        kernel(round(length(edges)/2)) = 1;
end

center = ceil(length(edges)/2);
for ii=1:size(spike_times_binned_further,2)
    s = conv(spike_times_binned_further(:,ii),kernel);
    s = s(center:end);
    s = s(1:length(spike_times_binned_further(:,ii)));
    smoothed_spike_times(:,ii) = s;
end

% Plot smoothed spike times
figure
for ii=1:size(smoothed_spike_times,2)
    plot(time_further_binned, ii-0.5-0.9*(smoothed_spike_times(:,ii)/...
            max(smoothed_spike_times(:,ii))))
    hold on
end
hold off
set(gca,'YDir','reverse')
yticks([0:1:(size(data,2)-1)])
yticklabels({'A';'B';'C'})
ylabel('Neuron')
ylim([min(ylim)-0.25 max(ylim)+0.25])
ylim_smoothed = ylim;

% Estimate CTM-DI
% First convert spike_times 1) from a cell to zero padded
% array and 2) from sec to ms
spike_times_ms_array = zeros(length(time),length(spike_times));
for ii=1:length(spike_times)
    spike_times_ms_array(1:length(spike_times{ii}),ii) = ...
        spike_times{ii}*(1e3);
end
[CMatrix, w_i, HMatrix] = CTM_spiketime_wrapper(spike_times_ms_array,M,bin_width,1);

% Compute + plot partial cross correlation
r_partial = nan(size(smoothed_spike_times,2),size(smoothed_spike_times,2),(2*M)+1);
connection_sign = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_sum = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_regular = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
connection_sign_regular_sum = zeros(size(smoothed_spike_times,2),size(smoothed_spike_times,2));
for ii=1:size(smoothed_spike_times,2)
    for jj=1:size(smoothed_spike_times,2)
        if ii~=jj
            [r_partial, r_regular] = partial_xcorr(smoothed_spike_times, ii, jj, M);
            r_partial_causal = r_partial(1:M);
            r_regular_causal = r_regular(1:M);
            
            I = find(max(abs(r_partial_causal))==abs(r_partial_causal));
            connection_sign(ii,jj) = sign(r_partial_causal(I));
            connection_sign_sum(ii,jj) = sum(r_partial_causal);
            
            I = find(max(abs(r_regular_causal))==abs(r_regular_causal));
            connection_sign_regular(ii,jj) = sign(r_regular_causal(I));
            connection_sign_regular_sum(ii,jj) = sum(r_regular_causal);            
        end
    end
end


% Estimate CCDI
%smoothed_spike_times(smoothed_spike_times>1)=1;
[DI, DI_list] = di_compute(smoothed_spike_times,M,C,boot_iter);

DI(DI<0) = 0;

DI_norm = DI*(1/log(2))./HMatrix;
DI_norm(DI_norm<0.01) = 0;
DI_thresholded = (DI_norm/(1/log(2))).*HMatrix;

thresh = eps;

DI_cond_post = di_compute_post(DI_thresholded,thresh,M,smoothed_spike_times,boot_iter);
DI_cond_post(DI_cond_post<0)=0;
DI_cond_post_norm = DI_cond_post*(1/log(2))./HMatrix;

DI_cond_post_norm(logical(eye(11))) = 0;
DI_norm(logical(eye(11))) = 0;

% SORTING
[DI_norm_SORTED, I_DI_norm_SORTED] = ...
    sort(DI_norm(:));
[DI_cond_post_norm_SORTED, I_DI_cond_post_norm_SORTED] = ...
    sort(DI_cond_post_norm(:));

save hh_11_neuron.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%