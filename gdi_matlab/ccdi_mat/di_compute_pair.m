function [DI, DI_list] = di_compute_pair(X,M,C,B,pairs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTES DI BETWEEN EACH CHANNEL (COLUMN) OF X USING WINDOWS (RANGES OF
% ROWS)
%
% INPUTS:
%   X - data with dim (observations)X(channels)
%   M - number of past samples to use for DI estimation
%   C - flag to indicate: (1) condition DI on rest of channels
%                         (0) don't condition DI on other nodes
%   B - number of bootstrap iterations
%   pairs - Px2 matrix, where each row is a pair and the first column is
%           the source nodes while the second column is the sink nodes. DI
%           will be computed from the source to the sink nodes. 
%           E.g., pairs = [1 2] means DI will be computed from 1 to 2.
% 
% OUTPUT: 
%   DI - matrix of DI values where the direction is DI from row to column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIMENSIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = size(X,1);              % Number of observations
    D = size(X,2);              % Number of channels
    W = floor(N/(M+1));         % Number of windows
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FORMAT FOR DI ESTIMATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X     = X(1:(W*(M+1)),:);     % Remove observations at end that can't use
    
    X_win = nan([W,M+1,D]);      % Window data indexed by window index
    for ii=1:W
        X_win(ii,:,:) = X((1:(M+1))+((ii-1)*(M+1)),:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE DI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    number_of_pairs = size(pairs, 1);               % # of pairs
    DI_vec = nan(number_of_pairs,1);         % Save DI as vec, unpack later
    DI_vec_list = nan(number_of_pairs,B);    % Similar to DI vec w/ all iter
    delete('ccdi_mat/progress/*')            % Clean progress folder
    for ii=1:number_of_pairs
        
        % Write file to progress folder to indicate program progress
        fileID = fopen(sprintf('ccdi_mat/progress/%i_of_%i.txt',...
                   ii,number_of_pairs),'w');
        fclose(fileID);
        
        % Select past of x and current values y
        x = squeeze(X_win(:,1:M,pairs(ii,1)));
        y = squeeze(X_win(:,M+1,pairs(ii,2)));
        
        % Check if conditioning on rest of channels or not
        if C==1
            z_array = X_win(:,1:M,:);
            z_array(:,:,pairs(ii,1)) = []; % Condition on all channels but x
            for jj=1:(D-1)
                z(:,(1:M)+((jj-1)*M)) = squeeze(z_array(:,:,jj));
            end
        elseif C==0  % Unconditional DI - use just past of y
            z = squeeze(X_win(:,1:M,pairs(ii,2)));
        else
            error('C is invalid value');
        end
            
        % Compute DI
        size(y)
        [DI_vec(ii), DI_vec_list(ii,:)] = ccdi_mat_matlab_fun(x,y,z,B);
    end
    
    % Repack MI into matrix
    DI = nan(D,D); % Initialize DI mat
    DI_list = nan(D,D,B); % Initialize DI_list mat
    for ii=1:number_of_pairs
        DI(pairs(ii,1), pairs(ii,2)) = DI_vec(ii);
        DI_list(pairs(ii,1), pairs(ii,2),:) = DI_vec_list(ii,:);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%