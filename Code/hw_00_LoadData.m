%% 0. Load the data

dt = 2;
year = 1845 : dt : 1903;
year = year';
m = length(year);

hare1 = [
    20
    20
    52
    83
    64
    68
    83
    12
    36
    150
    110
    60
    7
    10
    70
    ];
hare2 = [
    100
    92
    70
    10
    11
    137
    137
    18
    22
    52
    83
    18
    10
    9
    65
    ];
hare = [hare1; hare2];

lynx1 = [
    32
    50
    12
    10
    13
    36
    15
    12
    6
    6
    65
    70
    40
    9
    20
    ];
lynx2 = [
    34
    45
    40
    15
    15
    60
    80
    26
    18
    37
    50
    35
    12
    12
    25
    ];
lynx = [lynx1; lynx2];

clear hare1
clear hare2
clear lynx1
clear lynx2

%% Plotting

% Plot the data
% figure
% hold on
% plot(year, hare, 'b.--')
% plot(year, lynx, 'r.--')
% ax = gca;
% ax.FontSize = 12;
% yticks([0:25:150])
% xlabel('Year', 'fontsize', 20)
% ylabel('P', 'fontsize', 20)
% axis([1840 1910 0 150])
% ttl_str = sprintf('Population vs. Year');
% title(ttl_str, 'fontsize', 20)
% legend('Hare', 'Lynx', 'location', 'northeast')