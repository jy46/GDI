% INITIALIZE
clear all
close all

% PATHS
addpath ccdi_mat

% PARAMETERS
R  = 11;      % number of nodes
N  = (1e5)/5;    % number of samples
M  = 3;      % memory parameter for CCDI
p1 = 0.5; % parameters for connectivity
p2 = 0.75;
p3 = 0.5;
p4 = 0.5;
p5 = 0.75;
B = 2; % num bootstrap iter
 
% GENERATE
X = mvnrnd(zeros(1,R),eye(R),N);   % X has dim: NxR
X(2:end,1)  = sqrt(1-p1)*X(2:end,1)  + sqrt(p1)*X(1:(end-1),6);
X(2:end,3)  = sqrt(1-p1)*X(2:end,3)  + sqrt(p1)*X(1:(end-1),6);
X(2:end,5)  = sqrt(1-p1)*X(2:end,5)  - sqrt(p1)*X(1:(end-1),6);
X(2:end,4)  = sqrt(1-p2)*X(2:end,4)  + sqrt(p2)*X(1:(end-1),1);
X(2:end,2)  = sqrt(1-p3)*X(2:end,2)  ...%+ sqrt(p3/4)*X(1:(end-1),7)...
                                     - sqrt(p3/3)*X(1:(end-1),8)...
                                     + sqrt(p3/3)*X(1:(end-1),9)...
                                     + sqrt(p3/3)*X(1:(end-1),6);
X(2:end,11) = sqrt(1-p4)*X(2:end,11) + sqrt(p4)*X(1:(end-1),4);
X(2:end,10) = sqrt(1-p5)*X(2:end,10) + sqrt(p5/2)*X(1:(end-1),5)...
                                     + sqrt(p5/2)*X(1:(end-1),3);                               

% TRUE CONNECTIVITY
true_connectivity = zeros(R,R);
true_connectivity(6,[1 3 5]) = [p1 p1 -p1];
true_connectivity([6 9],2) = p3/3;
true_connectivity(8,2) = -p3/3;
true_connectivity([3 5],10) = p5/2;
true_connectivity(1,4) = p2;
true_connectivity(4,11) = p4;

% ESTIMATE SIGN                                 
[connection_sign, connection_sign_regular] = sign_inference(X,M)

% ESTIMATE DI
DI_uncond = di_compute(X,M,0,B);
DI_uncond(DI_uncond<0) = 0;
DI_uncond
DI_uncond_signed = connection_sign_regular.*DI_uncond
DI_uncond_signed(logical(eye(R)))=0

DI_cond = di_compute(X,M,1,B);
DI_cond(DI_cond<0) = 0;
DI_cond
DI_cond_signed = connection_sign.*DI_cond
DI_cond_signed(logical(eye(R)))=0

%% COMPUTE TRUE VALUES OF GDI
true_GDI = zeros(R,R);

true_GDI(6,1)  = 0.5*log(1+(p1/(1-p1)));
true_GDI(6,3)  = 0.5*log(1+(p1/(1-p1)));
true_GDI(6,5)  = -0.5*log(1+(p1/(1-p1)));

true_GDI(6,2)  =  0.5*log(1+((p3/3)/(1-p3)));
true_GDI(8,2)  = -0.5*log(1+((p3/3)/(1-p3)));
true_GDI(9,2)  =  0.5*log(1+((p3/3)/(1-p3)));

true_GDI(3,10) = 0.5*log(1+(((p5/2)*(1-p1))/(1-p5)));
true_GDI(5,10) = 0.5*log(1+(((p5/2)*(1-p1))/(1-p5)));

true_GDI(1,4)  = 0.5*log(1+(((p2)*(1-p1))/(1-p2)));
true_GDI(4,11) = 0.5*log(1+(((p4)*(1-p2))/(1-p4)));

%% PLOT RESULTS
figure

subplot(1,3,1)
imagesc(DI_uncond_signed)
colormap cool
colorbar
c_save(1,:) = caxis;
title('Estimated DI')

subplot(1,3,2)
imagesc(true_GDI)
colormap cool
colorbar
c_save(2,:) = caxis;
title('True GDI')

subplot(1,3,3)
imagesc(DI_cond_signed)
colormap cool
colorbar
c_save(3,:) = caxis;
title('Estimated GDI')

for ii=1:3
    subplot(1,3,ii)
    caxis([-max(abs(c_save(:))) max(abs(c_save(:)))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
