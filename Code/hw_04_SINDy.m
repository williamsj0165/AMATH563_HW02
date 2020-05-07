%% 4. SINDy

% Step 0. Build X and X_prm
X = [hare(1:end-1)'; lynx(1:end-1)']';
% % first-order forward finite difference for computing X_prm
X_prm = ([hare(2:end)'; lynx(2:end)'] - [hare(1:end-1)'; lynx(1:end-1)'] )' / dt;

% Step 1. Build library of candidate functions Theta_X from fncs of X
F = 22;
Theta_X = zeros((m-1),F);
for t = 1 :(m-1)
    Theta_X(t,1)  = 1;
    Theta_X(t,2)  = X(t,1);
    Theta_X(t,3)  = X(t,2);
    
    Theta_X(t,4)  = X(t,1).^2;
    Theta_X(t,5)  = X(t,1)*X(t,2);
    Theta_X(t,6)  = X(t,2).^2;
    
    Theta_X(t,7)  = X(t,1).^3;
    Theta_X(t,8)  = X(t,1).^2 * X(t,2);
    Theta_X(t,9)  = X(t,1) * X(t,2).^2;
    Theta_X(t,10)  = X(t,2).^3;
    
    Theta_X(t,11) = sin(X(t,1));
    Theta_X(t,12) = sin(X(t,2));
    Theta_X(t,13)= cos(X(t,1));
    Theta_X(t,14)= cos(X(t,2));
    
    Theta_X(t,15)= sin(2*X(t,1));
    Theta_X(t,16)= sin(2*X(t,2));
    Theta_X(t,17)= cos(2*X(t,1));
    Theta_X(t,18)= cos(2*X(t,2));
    
    Theta_X(t,19)= sin(3*X(t,1));
    Theta_X(t,20)= sin(3*X(t,2));
    Theta_X(t,21)= cos(3*X(t,1));
    Theta_X(t,22)= cos(3*X(t,2));
end

% Step. 2 Compute Xi matrix
n = 2;
lambda = 2.5;
Xi_stl = sparsifyDynamics(Theta_X,X_prm,lambda,n);


% Step 3. Build symbolic Theta_x
Theta_x = zeros(t,F);
x_dot_stl = zeros(2,(m-1));
for t = 1 : (m-1)
    for k = 1 : 2
        Theta_x(t,:) = [ 1 X(t,1) X(t,2) X(t,1)^2 X(t,1)*X(t,2) X(t,2)^2 X(t,1)^3 X(t,1)^2*X(t,2) X(t,1)*X(t,2)^3 X(t,2)^3 sin(X(t,1)) cos(X(t,1)) sin(X(t,2)) cos(X(t,2)) sin(2*X(t,1)) cos(2*X(t,1)) sin(2*X(t,1)) sin(2*X(t,2)) sin(3*X(t,1)) cos(3*X(t,1)) sin(3*X(t,2)) cos(3*X(t,2))];
        x_dot_stl(k,t) = Theta_x(t,:) * Xi_stl(:,k);
    end
end

t_fut = 50;
x_pred_stl = zeros(t_fut,2);
x_pred_stl(1:m-1,:) = X;
for t = m : t_fut
    x_pred_stl(t,1) = SINDy_lib( x_pred_stl(t-1,:), Xi_stl(:,1)) * dt + x_pred_stl(t-1,1);
    x_pred_stl(t,2) = SINDy_lib( x_pred_stl(t-1,:), Xi_stl(:,2)) * dt + x_pred_stl(t-1,2);
end

hare_zro = 0;
lynx_zro = 0;
for f = 1 : F
    if abs(Xi_stl(f,1)) < 1E-6
        hare_zro = hare_zro + 1;
    end
    if abs(Xi_stl(f,2)) < 1E-6
        lynx_zro = lynx_zro + 1;
    end
end

K_1 = F - hare_zro;
K_2 = F - lynx_zro;


%% Plotting

figure
hold on
plot(Xi_stl(:,1),'bo--')
plot(Xi_stl(:,2),'ro--')
xticks([0:1:F])
xlabel('Term', 'fontsize', 20)
ylabel('\xi', 'fontsize', 20)
axis([0 23 -15 30])
ttl_str = sprintf('Weights of Library Functions\n\\lambda = %0.4f SINDy Prediction\nK_1 = %d, K_2 = %d',lambda,K_1,K_2);
title(ttl_str, 'fontsize', 20)
legend('Hare', 'Lynx', 'location', 'northeast')

[[1:F]' Xi_stl]


%% Plotting

% Plot data against SINDy predictions for HARE
figure
hold on
plot(year, hare, 'b*--')
plot(1845:2:1845+2*(t_fut-1), x_pred_stl(:,1), 'c.-')
ax = gca;
ax.FontSize = 12;
yticks([0:25:175])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 175])
ttl_str = sprintf('Hare Population vs. Year\n\\lambda = %0.4f SINDy Prediction to 1943',lambda);
title(ttl_str, 'fontsize', 20)
legend('Data', 'SINDy', 'location', 'northeast')

% Plot data against SINDy predictions for LYNX
figure
hold on
plot(year, lynx, 'r*--')
plot(1845:2:1845+2*(t_fut-1), x_pred_stl(:,2), 'm.-')
ax = gca;
ax.FontSize = 12;
yticks([0:25:100])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1950 0 100])
ttl_str = sprintf('Lynx Population vs. Year\n\\lambda = %0.4f SINDy Prediction to 1943',lambda);
title(ttl_str, 'fontsize', 20)
legend('Data', 'SINDy', 'location', 'northeast')



