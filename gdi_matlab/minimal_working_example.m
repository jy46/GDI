addpath ccdi_mat

num_samples = 5000;
mean_vec = [0,0,0];
cov_mat  = eye(3);

X = mvnrnd(mean_vec, cov_mat, num_samples);

X(3:end,2) = X(3:end,2) + X(1:(end-2),1);
X(3:end,3) = X(3:end,3) + X(1:(end-2),2);

M = 4;
B = 10;

X_DI = DI(X,M,B)
X_GDI = GDI(X,M,B)
