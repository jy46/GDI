function [connection_sign, connection_sign_regular] = sign_inference(X,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM SIGN INFERENCE OF RELATIONSHIPS BETWEEN COLUMNS OF X
%
% INPUT:
%   X: Input data with dim (sample)x(channel)
%   M: History length, i.e. number of past samples to use
% OUTPUT:
%   connection_sign: Sign inferred for relationships between columns of X
%                    using partial correlations. connection_sign has shape 
%                    (channel)x(channel) and is the sign of the relationship
%                    from the row channel to the column channel.
%   connection_sign_regular: Sign inferred for relationships between columns of X
%                    using regular correlations. connection_sign_regular has shape 
%                    (channel)x(channel) and is the sign of the relationship
%                    from the row channel to the column channel.
%
% Copyright (C) 2020 Joseph Young - see GPLv2_note.txt for full notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_partial = nan(size(X,2),size(X,2),(2*M)+1);
    connection_sign = zeros(size(X,2),size(X,2));
    r_regular = nan(size(X,2),size(X,2),(2*M)+1);
    connection_sign_regular = zeros(size(X,2),size(X,2));
    for ii=1:size(X,2)
        for jj=1:size(X,2)
            if ii~=jj
                [r_partial, r_regular] = partial_xcorr(X, ii, jj, M);

                r_partial_causal = r_partial(1:M);
                r_regular_causal = r_regular(1:M);

                I = find(max(abs(r_partial_causal))==abs(r_partial_causal));
                connection_sign(ii,jj) = sign(r_partial_causal(I));

                I = find(max(abs(r_regular_causal))==abs(r_regular_causal));
                connection_sign_regular(ii,jj) = sign(r_regular_causal(I));                 
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
