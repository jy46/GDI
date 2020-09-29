%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INCREASING SAMPLE SIZE FIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE
clear all
close all


% FLDRS
addpath ccdi_mat


% PARAMETERS
sample_sizes = [100 500 1000 2000 5000 1e4 2e4 1e5 2e5]; % x axis
%sample_sizes = [2e3 5e3];
dim_sizes    = 2+[0 10 20 50 100]; % multiple lines, each for fixed dim
%dim_sizes = 2+[0 10];
boot_iter    = 10; % boot iterations per computation
num_runs     = 50; % number of runs for each data point
sim_type     = 'discrete'; % 'discrete' or 'continuous'
    rho = 0.6; % covariance for continuous

    p  = 0.1; % for discrete, effectiveness of synaptic transmission
    Px = 0.3; % prob of x spiking in a given time bin
    

% LOOP THROUGH SIMULATIONS
% https://www.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-
% all-possible-combinations-of-given-vectors-in-matlab#answer_252633
elements = {dim_sizes, sample_sizes, 1:num_runs}; %cell array with N vectors to combine
combinations = cell(1, numel(elements)); %set up the varargout result
[combinations{:}] = ndgrid(elements{:});
combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); 
iter_mat = [combinations{:}]; % NumberOfCombinations by N matrix. Each row unique.

DI = nan(size(iter_mat,1),1);
DI_list = nan(size(iter_mat,1),boot_iter);
delete('ccdi_mat/pair_progress/*')            % Clean progress folder
parfor ii=1:size(iter_mat,1)
    
    % Write file to progress folder to indicate program progress
    fileID = fopen(sprintf('ccdi_mat/pair_progress/%i_of_%i.txt',...
               ii,size(iter_mat,1)),'w');
    fclose(fileID);
    
    % Current parameters
    dim_size    = iter_mat(ii,1);
    sample_size = iter_mat(ii,2);
    run         = iter_mat(ii,3);
    
    if strcmp(sim_type,'continuous')
        
        % COVARIANCE
        C = eye(dim_size);
        C(1,2) = rho;   % correlation between dim 1 and dim 2
        C(2,1) = rho;

        % GENERATE DATA
        X = mvnrnd(zeros(1,dim_size),...  % mean
                   C,...                             % covariance
                   sample_size);       % num samples
        
    elseif strcmp(sim_type,'discrete')
        
        X = randi(2,sample_size,dim_size)-1;
        X(:,1) = binornd(1,Px,sample_size,1);
        X(:,2) = bsc(X(:,1),p);
        
    end
    
    % SHIFT FOR CAUSALITY: DI IS FROM 1 TO 2; OTHER DIMS RANDOM
    X(2:end,2) = X(1:(end-1),2);
    
    % COMPUTE DI
    [DI_raw, DI_list_raw] = ...
        di_compute_pair(X,1,1,boot_iter,[1 2]);
    DI(ii) = DI_raw(1,2);
    DI_list(ii,:) = DI_list_raw(1,2,:);
             
end


% UNROLL AND PUT INTO DI MAT
DI_mat = nan(length(dim_sizes),length(sample_sizes),num_runs);
DI_list_mat = nan(length(dim_sizes),length(sample_sizes),num_runs,boot_iter);
for ii=1:size(iter_mat,1)
    
    dim_size    = iter_mat(ii,1);
    sample_size = iter_mat(ii,2);
    run         = iter_mat(ii,3);
    
    DI_mat(dim_sizes==dim_size,sample_sizes==sample_size,run) = DI(ii);
    DI_list_mat(dim_sizes==dim_size,sample_sizes==sample_size,run,:) = DI_list(ii,:);
end


% TRUE DI
if strcmp(sim_type,'continuous')
    true_DI = -0.5*log(1-(rho^2));
elseif strcmp(sim_type,'discrete')
    X = binornd(1,Px,max(sample_sizes),1);
    Y = bsc(X,p);
    Y(2:end) = Y(1:end-1);
    py1 = sum(Y)/max(sample_sizes);
    true_DI = (-py1*log(py1)-(1-py1)*log(1-py1))-(-p*log(p)-(1-p)*log(1-p));
end

%% SAVE
save(sprintf('scaling_figs/%s_saved.mat',sim_type))

%% PLOT AGAINST SAMPLE SIZE
figure
subplot(1,3,1)
for ii=1:size(DI_mat,1)
    errorbar(sample_sizes/2,mean(DI_mat(ii,:,:),3),var(DI_mat(ii,:,:),[],3))
    set(gca,'XScale','log')
    %errorbarlogx()
    hold on
end
line([min(xlim) max(xlim)],[true_DI true_DI],'Color','k','LineStyle','--')
hold off
legend([arrayfun(@(x) sprintf('d_z = %i',x-2),dim_sizes,'UniformOutput',false),...
       {'True DI'}],'Location','southeast')
xlabel('N')
ylabel('DI')
title('Estimated DI vs. Sample Size')
box off
ylim([0 max(ylim)])


% PLOT AGAINST DIM
subplot(1,3,2)
for ii=1:size(DI_mat,2)
    errorbar(dim_sizes-2,mean(DI_mat(:,ii,:),3),var(DI_mat(:,ii,:),[],3))
    hold on
end
line([min(xlim) max(xlim)],[true_DI true_DI],'Color','k','LineStyle','--')
hold off
legend([arrayfun(@(x) sprintf('N = %i',x),sample_sizes/2,'UniformOutput',false),...
       {'True DI'}],'Location','southeast')
xlabel('d_z')
ylabel('DI')
title('Estimated DI vs. Dimensionality of Z')
box off
ylim([0 max(ylim)])


% PLOT AGAINST BOOT
subplot(1,3,3)
for ii=1:min(num_runs,5)
    plot(1:boot_iter,cumsum(squeeze(DI_list_mat(1,5,ii,:)))'./(1:boot_iter))
    hold on
end
line([min(xlim) max(xlim)],[true_DI true_DI],'Color','k','LineStyle','--')
hold off
legend([arrayfun(@(x) sprintf('Run %i',x),1:min(num_runs,5),'UniformOutput',false),...
       {'True DI'}],'Location','southeast')
xlabel('# Bootstrap Iterations')
ylabel('DI')
title('Estimated DI vs. Bootstrap Iterations')
box off
ylim([0 max(ylim)])


%% SAVE





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%