function [r_partial, r] = partial_xcorr(X_raw, reference, target, max_tau)
% Partial cross correlation from reference (index of channel) to target 
% (index of other channel). 
% X_raw is input data in shape (observations)x(channels)
% max_tau is the maximum integer to shift positively and negatively for the
% partial cross correlation. All other channels are conditioned on in this
% analysis.

    % Initialize
    tau_vec = (-max_tau:max_tau);
    conditional_regressor_columns = 1:size(X_raw,2);
    conditional_regressor_columns([reference target]) = [];
    
    Z = X_raw(:,conditional_regressor_columns);
    X = X_raw(:,reference);
    Y = X_raw(:,target);
    
    r_partial = nan(size(tau_vec));
    r         = nan(size(tau_vec));
    
    % Loop through lags
    for ii=1:length(tau_vec)
        
        % Reset so zeros will pad things
        Y_lag = zeros(size(Y));
        
        % Retrieve lagged values
        indices = (1:length(Y)) - tau_vec(ii);
        if tau_vec(ii)==0
            Y_lag = Y;
        elseif tau_vec(ii)<0
            indices(indices>length(Y)) = [];
            Y_lag(1:length(indices)) = Y(indices);
        elseif tau_vec(ii)>0
            indices(indices<1) = [];
            Y_lag(end-length(indices)+1:end) = Y(indices);
        else
            error('ERROR');
        end
        
        % Compute partial correlation at lag
        r_partial(ii) = partialcorr(X,Y_lag,Z);
        
        % Also compute correlation at lag
        r_mat = corr([X Y_lag]);
        r(ii) = r_mat(1,2);
    end

end