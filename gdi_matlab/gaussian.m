% INITIALIZE
clear all
close all

% PATHS
addpath ccdi_mat

% PARAMETERS
R  = 11;      % number of nodes
N  = 1e5;    % number of samples
M  = 3;      % memory parameter for CCDI
p1 = 0.5;
p2 = 0.75;
p3 = 0.5;
p4 = 0.5;
p5 = 0.75;
B = 100; % num bootstrap iter
 
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

% ESTIMATE DI
DI_uncond = di_compute(X,M,0,B);
DI_uncond(DI_uncond<0) = 0;
DI_uncond
DI_uncond_signed = connection_sign_regular.*DI_uncond
DI_uncond_signed(abs(DI_uncond_signed)<0.05)=0
DI_uncond_signed(logical(eye(R)))=0

DI_cond = di_compute(X,M,1,B);
DI_cond(DI_cond<0) = 0;
DI_cond
DI_cond_signed = connection_sign.*DI_cond
DI_cond_signed(abs(DI_cond_signed)<0.05)=0
DI_cond_signed(logical(eye(R)))=0

save gaussian_example.mat

% CREATE GRAPHS
% COND
[s,t] = find(DI_cond_signed~=0);
c_weights = DI_cond_signed(sub2ind(size(DI_cond_signed), s, t));
G_DI_cond_signed = digraph(s,t,c_weights);

% UNCOND
[s,t] = find(DI_uncond_signed~=0);
u_weights = DI_uncond_signed(sub2ind(size(DI_uncond_signed), s, t));
G_DI_uncond_signed = digraph(s,t,u_weights);

%%
close all

scale_factor = 10;

figure
h_c = plot(G_DI_cond_signed,'LineWidth',scale_factor*abs(G_DI_cond_signed.Edges.Weight));
box off
title('Conditioned')

xlim([0 3])
ylim([0 3])

h_c.NodeColor = [141,211,199]/255;
h_c.MarkerSize = 30;
h_c.EdgeColor = 'k';
h_c.EdgeAlpha = 1;
h_c.ArrowPosition = 0.6;
h_c.NodeLabel = [];
set(gca,'XColor','none')
set(gca,'YColor','none')
axis square
h_c.ArrowSize = 15;


h_c.XData(6) = mean(xlim);
h_c.YData(6) = mean(ylim);

h_c.XData(2) = mean(xlim) + 0.2*diff(xlim);
h_c.YData(2) = h_c.YData(6);

h_c.XData(7) = h_c.XData(2);
h_c.YData(7) = h_c.YData(6) + 0.2*diff(ylim);

h_c.XData(8) = h_c.XData(2) + 0.2*diff(xlim);
h_c.YData(8) = h_c.YData(2);

h_c.XData(9) = h_c.XData(2);
h_c.YData(9) = h_c.YData(6) - 0.2*diff(ylim);

h_c.XData(5) = h_c.XData(6);
h_c.YData(5) = h_c.YData(6) - 0.2*diff(ylim);

[rotpoint] = rotate_in_2D([h_c.XData(5)-mean(xlim) h_c.YData(5)-mean(ylim)]',-5);
h_c.XData(5) = rotpoint(1)+mean(xlim);
h_c.YData(5) = rotpoint(2)+mean(ylim);

[rotpoint] = rotate_in_2D([h_c.XData(5)-mean(xlim) h_c.YData(5)-mean(ylim)]',-45);
h_c.XData(3) = rotpoint(1)+mean(xlim);
h_c.YData(3) = rotpoint(2)+mean(ylim);

h_c.XData(10) = (h_c.XData(3) + h_c.XData(5))-mean(xlim);
h_c.YData(10) = (h_c.YData(3) + h_c.YData(5))-mean(ylim);

[rotpoint] = rotate_in_2D([h_c.XData(5)-mean(xlim) h_c.YData(5)-mean(ylim)]',-170);
h_c.XData(1) = 1.1*rotpoint(1)+mean(xlim);
h_c.YData(1) = 1.1*rotpoint(2)+mean(ylim);

[rotpoint] = rotate_in_2D([h_c.XData(1)-mean(xlim) h_c.YData(1)-mean(ylim)]',45);
h_c.XData(4) = rotpoint(1)+h_c.XData(1);
h_c.YData(4) = rotpoint(2)+h_c.YData(1);

[rotpoint] = rotate_in_2D([h_c.XData(1)-mean(xlim) h_c.YData(1)-mean(ylim)]',90);
h_c.XData(11) = rotpoint(1)+h_c.XData(4);
h_c.YData(11) = rotpoint(2)+h_c.YData(4);


for ii=1:11
    if ii<10
        text(h_c.XData(ii)-0.015*diff(xlim),h_c.YData(ii),num2str(ii),...
             'FontSize', 12)
    else
        text(h_c.XData(ii)-0.03*diff(xlim),h_c.YData(ii),num2str(ii),...
             'FontSize', 12)
    end
end

% viscircles([h_c.XData' h_c.YData'],0.19*ones(11,1),'Color','k',...
%            'EnhanceVisibility',0)



figure
h_u = plot(G_DI_uncond_signed,'LineWidth',scale_factor*abs(G_DI_uncond_signed.Edges.Weight));
box off

xlim([0 3])
ylim([0 3])

h_u.XData = h_c.XData;
h_u.YData = h_c.YData;
title('Pairwise')
h_u.NodeColor = [141,211,199]/255;
h_u.MarkerSize = 30;
h_u.EdgeColor = 'k';
h_u.EdgeAlpha = 1;
h_u.ArrowPosition = 0.6;
h_u.NodeLabel = [];
set(gca,'XColor','none')
set(gca,'YColor','none')
axis square
h_u.ArrowSize = 15;


for ii=1:11
    if ii<10
        text(h_c.XData(ii)-0.015*diff(xlim),h_c.YData(ii),num2str(ii),...
             'FontSize', 12)
    else
        text(h_c.XData(ii)-0.03*diff(xlim),h_c.YData(ii),num2str(ii),...
             'FontSize', 12)
    end
end


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








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
