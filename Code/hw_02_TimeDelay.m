%% 2. Time-Delay Embedded DMD Model

% Step 0. Build X, X_prm, and the Hankel mtx H
m=length(year);

D = 25; % degree of time-delay embedding
W = 5; % width of time-delay embedding
if W + D > 30
    fprintf('\nError! W + D >= 30\n')
    fprintf('\nAborting.\n\n')
    return
end
H = zeros(D*2,W);
H_prm = zeros(D*2,W);
t=1;

for d = 1 : 2 : D*2
        H(d,:)   = hare(t:W+(t-1));
        H(d+1,:) = lynx(t:W+(t-1));
        
        H_prm(d,:)   = hare(t+1:W+t);
        H_prm(d+1,:) = lynx(t+1:W+t);
        
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

%% Plotting

% Plot data against DMD predictions for HARE
figure
hold on
plot(year, hare, 'b.--')
plot(1845:2:1845+2*(t_fut-1), h(1,:), 'c*--')
ax = gca;
ax.FontSize = 12;
yticks([0:25:150])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 150])
ttl_str = sprintf('Hare Population vs. Year\n(m_c,m_o) = (%0.d,%0.d) Time-Delay DMD Prediction to 1943',W,D);
title(ttl_str, 'fontsize', 20)
legend('Data', 'Time-Delay DMD', 'location', 'northeast')

% Plot data against DMD predictions for LYNX
figure
hold on
plot(year, lynx, 'r.--')
plot(1845:2:1845+2*(t_fut-1), h(2,:), 'm*--')
ax = gca;
ax.FontSize = 12;
yticks([0:25:100])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 100])
ttl_str = sprintf('Lynx Population vs. Year\n(m_c,m_0) = (%0.d,%0.d) Time-Delay DMD Prediction to 1943',W,D);
title(ttl_str, 'fontsize', 20)
legend('Data', 'Time-Delay DMD', 'location', 'northeast')

