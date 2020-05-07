%% 3. Lotka-Volterra Model

t1_tm = [1845 1943];
y0 = [20; 32];   
[t,y] = ode23(@LotkaVolterra,t1_tm,y0);

% Guessed values of parameters are:
b = 0.5;
p = 0.04;
r = 0.02;
d = 1;


%% A form for KL divergence
dt_LV = t(2) - t(1);

t_KL = 1845 : 2 : 1903;
y_KL = zeros(2,m);
for tstep = 1 : length(t_KL)
    [M,I] = min(abs(t-t_KL(tstep)));
    y_KL(:,tstep) = y(I,:);
end


%% Plotting

% Plot data against Lotka-Volterra predictions for HARE
figure
hold on
plot(year, hare, 'b.--')
plot(t, y(:,1), 'c')
% plot(t_KL,y_KL(2,:), '*k--', 'linewidth', 1)
ax = gca;
ax.FontSize = 12;
yticks([0:25:150])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1943 0 150])
ttl_str = sprintf('Hare Population vs. Year\nLotka-Volterra Prediction to 1943\nb = %0.2f, p = %0.2f, r = %0.2f, d = %0.2f',b,p,r,d);
title(ttl_str, 'fontsize', 20)
legend('Data', 'Lotka-Volterra', 'location', 'northeast')

% Plot data against Lotka-Volterra predictions for LYNX
figure
hold on
plot(year, lynx, 'r.--')
plot(t, y(:,2), 'm')
% plot(t_KL,y_KL(2,:), '*k--', 'linewidth', 1)
ax = gca;
ax.FontSize = 12;
yticks([0:25:150])
xlabel('Year', 'fontsize', 20)
ylabel('P', 'fontsize', 20)
axis([1840 1943 0 100])
ttl_str = sprintf('Lynx Population vs. Year\nLotka-Volterra Prediction to 1943\nb = %0.2f, p = %0.2f, r = %0.2f, d = %0.2f',b,p,r,d);
title(ttl_str, 'fontsize', 20)
legend('Data', 'Lotka-Volterra', 'location', 'northeast')











