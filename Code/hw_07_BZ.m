%% 7. BZ Data

% loads the variable BZ_tensor
load BZ.mat
[m,n,k]=size(BZ_tensor);

fprintf('\nFinished loading data.\n\n')


%% Sample the pixel data

d_m_new = 50;
d_n_new = 50;
m_new = 1 : d_m_new : 351;
n_new = 1 : d_n_new : 451;
M = length(m_new)
N = length(n_new)

fprintf('\nFinished sampling data.\n\n')

%% Code provided in handout for viewing data
% Apparently used for viewing images as as video

[m,n,k]=size(BZ_tensor); % x vs y vs time data
for j=1:k
    A=BZ_tensor(:,:,j);
    pcolor(A), shading interp, pause(0.2)
end


%% Rescale BZ_tensor


max_BZ = max(max(max(BZ_tensor)));
min_BZ = min(min(min(BZ_tensor)));

BZ_tensor = (BZ_tensor - min_BZ) * 1/max_BZ;
% Checks
% max(max(max(BZ_tensor)))
% min(min(min(BZ_tensor)))

fprintf('\nFinished rescaling BZ_data.\n\n')


%% Form X matrices from BZ Data

% Step 0.
X_BZ = zeros(M*N,k-1);
for j = 1 : k-1
    fprintf('j_bz  = %d\n',j)
    X_BZ(:,j) = reshape(BZ_tensor(1:d_m_new:m,1:d_n_new:n,j),M*N,1);
end

X_BZ_prm = zeros(M*N,k-1);
for j = 2 : k
    fprintf('j_bz_prm  = %d\n',j)
    X_BZ_prm(:,j-1) = reshape(BZ_tensor(1:d_m_new:m,1:d_n_new:n,j),M*N,1);
end

fprintf('\nFinished forming X_BZ, X_BZ_prm.\n\n')


%% Initial Steps of KL Divergences

X = [X_BZ X_BZ_prm(:,end)];

% Does it matter if these are different?
range_BZ = 0 : 0.01 : 1;
l_BZ = length(range_BZ);

f_BZ = zeros(M*N,l_BZ);
for p = 1 : M*N
    fprintf('p_f_bz  = %d\n',p)
    f_BZ(p,:) = hist(X(p,:),range_BZ)+0.01;
    f_BZ(p,:) = f_BZ(p,:) / trapz(range_BZ,f_BZ(p,:));
end

fprintf('Finished with initial KL divergence.\n\n')


%% DMD

% Step 1. SVD of X
[U_BZ,Sigma_BZ,V_BZ] = svd(X_BZ,'econ');

% Step 2. Compute reduced form of A
A_tld_BZ = U_BZ' * X_BZ_prm * V_BZ * inv(Sigma_BZ);

% Step 3. Eigendecomposition of A_tld
[Q_BZ,Lambda_BZ] = eig(A_tld_BZ);

% Step 4. DMD Modes
Phi_BZ = X_BZ_prm * V_BZ * inv(Sigma_BZ) * Q_BZ;

% Step 5. DMD Expansion
b_BZ = Phi_BZ \ X_BZ(:,1);

% Step 6. Predict Future States
t_fut_BZ = k;
x_BZ = zeros(M*N,t_fut_BZ);
for tstep = 1 : t_fut_BZ
    fprintf('tstep_dmd_bz  = %d\n',tstep)
    x_BZ(:,tstep) = Phi_BZ * Lambda_BZ^(tstep - 1) * b_BZ;
end

fprintf('Finished with DMD_BZ.\n\n')


%% DMD KL Divergence, AIC & BIC Scores
[m,n,k]=size(BZ_tensor);

% KL Divergence
g1_BZ = zeros(M*N,l_BZ);
Int1_BZ = zeros(M*N,l_BZ);
I1_BZ = zeros(M*N,1);
for p = 1 : M*N
    fprintf('p_DMD_KL  = %d\n',p)
    g1_BZ(p,:) = hist(x_BZ(p,:),range_BZ)+0.01;
    g1_BZ(p,:) = g1_BZ(p,:) / trapz(range_BZ,g1_BZ(p,:));
    Int1_BZ(p,:) = f_BZ(p,:).*log(f_BZ(p,:)./g1_BZ(p,:));
    I1_BZ(j) = trapz(range_BZ,Int1_BZ(p,:));
end

I1_BZ_avg = mean(I1_BZ)

fprintf('Finished with KL divergence.\n\n')


% AIC and BIC Scores

%% Time-Delay Embedding
T = k;

D = 800; % degree of time-delay embedding
W = 400; % width of time-delay embedding
if W + D > T
    fprintf('\nError! W + D >= %d\n',T)
    fprintf('\nAborting.\n\n')
    return
end
H = zeros(D*2,W);
H_prm = zeros(D*2,W);
t=1;

for d = 1 : 2 : D*2
        for var = 1 : M*N
            H(d,var)   = X(var,t:W+(t-1));
            H_prm(d,var)   = X(var,t+1:W+t);
        end
        t=t+1;
end

% Step 1. SVD of H
[U_H,Sigma_H,V_H] = svd(H,'econ');

% Step 2. Compute reduced form of A
A_tld_H = U_H' * H_prm * V_H * inv(Sigma_H);

% Step 3. Eigendecomposition of A_tld
[Q_H,Lambda_H] = eig(A_tld_H);

% Step 4. DMD Modes
Phi_H = H_prm * V_H * inv(Sigma_H) * Q_H;

% Step 5. DMD Expansion
b_H = Phi_H \ H(:,1);

% Step 6. Predict Future States
t_fut = 50;
h = zeros(D*2,t_fut);
for k = 1 : t_fut
    h(:,k) = Phi_H * Lambda_H^(k - 1) * b_H;
end

h = real(h);


%% TDE KL Divergence









