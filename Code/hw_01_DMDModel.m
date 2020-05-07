%% 1. DMD Model

% Step 0. Build X and X_prm
X = [hare(1:end-1)'; lynx(1:end-1)'];
X_prm = [hare(2:end)'; lynx(2:end)'];

% Step 1. SVD of X
[U,Sigma,V] = svd(X,'econ');

% Step 2. Compute reduced form of A
A_tld = U' * X_prm * V * inv(Sigma);

% Step 3. Eigendecomposition of A_tld
[Q,Lambda] = eig(A_tld);

% Step 4. DMD Modes
Phi = X_prm * V * inv(Sigma) * Q;

% Step 5. DMD Expansion
b = Phi \ X(:,1);

% Step 6. Predict Future States
t_fut = 50;
x = zeros(2,t_fut);
for k = 1 : t_fut
    x(:,k) = Phi * Lambda^(k - 1) * b;
end

% % Step 7. Convert to continuous time
% Do I need to do this?
% I don't think so.
% dt_local = 2;
% Omega = log(Lambda) / dt_local;
% Omega(2,1) = -Inf;
% Omega(1,2) = -Inf;
% t = 1845 : dt_local : 1943;
% t = t - 1845;
% x_alt = zeros(2,length(t));
% for tstep = 1 : length(t)
%     x_alt(:,tstep) = Phi * exp(Omega*t(tstep)) * b;
% end
% t = t +1845;

%% Plotting

% Plot data against DMD predictions for HARE
figure
hold on
plot(year, hare, 'b.--')
plot(1845:2:1845+2*(t_fut-1), x(1,:), 'c*--')
% plot(t, x_alt(1,:), 'g.--')
ax = gca;
ax.FontSize = 12;
yticks([0:25:150])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 150])
ttl_str = sprintf('Hare Population vs. Year\nDMD Prediction to 1943');
title(ttl_str, 'fontsize', 20)
legend('Data', 'DMD', 'location', 'northeast')
% legend('Data', 'DMD', ' 7.30', 'location', 'northeast')

% Plot data against DMD predictions for LYNX
figure
hold on
plot(year, lynx, 'r.--')
plot(1845:2:1845+2*(t_fut-1), x(2,:), 'm*--')
% plot(t, x_alt(2,:), 'g.--')
ax = gca;
ax.FontSize = 12;
yticks([0:25:100])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 100])
ttl_str = sprintf('Lynx Population vs. Year\nDMD Prediction to 1943');
title(ttl_str, 'fontsize', 20)
legend('Data', 'DMD', 'location', 'northeast')
% legend('Data', 'DMD', ' 7.30', 'location', 'northeast')