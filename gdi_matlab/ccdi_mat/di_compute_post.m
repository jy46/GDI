function [DI_cond] = di_compute_post(DI_uncond,thresh,M,X,boot_iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDS INPUT UNCONDITIONED DI MATRIX BY THRESH, THEN COMPUTES
% CONDITIONAL DI BASED ON ABOVE THRESHOLD VALUES IN EACH COLUMN.
%
% INPUTS
%   DI_uncond: DxD matrix of unconditional DI.
%   thresh   : Threshold to apply to DI_uncond. DI_cond will only be based
%              on conditioning on channels with DI_uncond values >=thresh.
%   M        : Memory for DI.
%
% OUTPUTS
%   DI_cond  : DxD matrix of conditional DI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % INITIALIZE
    D       = size(DI_uncond,2); % Number of processes/channels
    DI_cond = DI_uncond; % DI_cond will be update in the following loop
    
    % COMPUTE CONDITIONAL DI
    for col=1:D
        
        % PRINT CURRENT COL
        fprintf('\n\n\n\n\n\n\n\n\n\n\n\n')
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('---------%d ---------\n',col)
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('---------------------\n')
        fprintf('\n\n\n\n\n\n\n\n\n\n\n\n')
        
        % CURRENT COLUMN VALUES
        current_col_values = DI_uncond(:,col);
        
        % FIND CHANNELS IN COLUMN WITH DI_UNCOND>=THRESH
        chan_with_DI_above_thresh = (find(current_col_values>=thresh))';
        
        % REMOVE CHANNEL FROM LIST IF CURRENT CHANNEL (col)
        chan_with_DI_above_thresh(chan_with_DI_above_thresh==col)=[];
        
        % COMPUTE DI_COND FOR THESE CHANNELS, OR LEAVE ALONE IF NO
        % CONDITIONING IS TO BE DONE (I.E. ONLY ONE CONNECTION)
        if length(chan_with_DI_above_thresh)>1
            [DI_col,~] = di_compute(X(:,[chan_with_DI_above_thresh col]),...
                                    M,1,boot_iter);
            DI_cond(chan_with_DI_above_thresh,col) = DI_col(1:(end-1),end);
        end
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
