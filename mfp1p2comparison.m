clear
close all
clc

rng(5); % Set random seed for reproducibility
p1 = 1;
p2 = 0.75;
sz = 2048;
k = 10;
q = -10:0.1:10;

fontsz = 24; % Global font size
%% Square MFp1p2 Spectrum Comparison

[initmat, alpha_theory, f_theory] = mfp1p2gen(sz, p1, p2, k);
[Dq_rec,myalpha_rec,falpha_rec] = mfrectanglebinarized(initmat, 0, q, 0);
[Dq_radial,myalpha_radial,falpha_radial] = mfthetacoordinate(initmat, 0, q, 0);
[Dq_polar,myalpha_polar,falpha_polar] = mfradialellipse(initmat, 0, q, 0);

% Plotting
h_fig_plot = figure(1); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', [0.85 0.325 0.098], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Rectangular'); % A warm orange-red
plot(myalpha_radial, falpha_radial, 'Color', [0.494 0.184 0.556], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Radial');   % A deep purple
plot(myalpha_polar, falpha_polar, 'Color', [0.301 0.745 0.933], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Polar');     % A vibrant light blue

hold off; 

xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', fontsz);
legend('Location', 'bestoutside','FontSize',16);
ax = gca;
ax.FontSize = 16; 

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 700, 500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font
%% Radial MFp1p2
[initmat, alpha_theory, f_theory] = mfp1p2radial(sz, p1, p2, k);
[Dq_rec,myalpha_rec,falpha_rec] = mfrectanglebinarized(initmat, 0, q, 0);
[Dq_radial,myalpha_radial,falpha_radial] = mfthetacoordinate(initmat, 0, q, 0);
[Dq_polar,myalpha_polar,falpha_polar] = mfradialellipse(initmat, 0, q, 0);

% Plotting
h_fig_plot = figure(2); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', [0.85 0.325 0.098], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Rectangular'); % A warm orange-red
plot(myalpha_radial, falpha_radial, 'Color', [0.494 0.184 0.556], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Radial');   % A deep purple
plot(myalpha_polar, falpha_polar, 'Color', [0.301 0.745 0.933], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Polar');     % A vibrant light blue

hold off; 

xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', fontsz);
legend('Location', 'bestoutside','FontSize',16);
ax = gca;
ax.FontSize = 16; 

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 700, 500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font

%% Polar MFp1p2
[initmat, alpha_theory, f_theory] = mfp1p2polar(sz, p1, p2, k);
[Dq_rec,myalpha_rec,falpha_rec] = mfrectanglebinarized(initmat, 0, q, 0);
[Dq_radial,myalpha_radial,falpha_radial] = mfthetacoordinate(initmat, 0, q, 0);
[Dq_polar,myalpha_polar,falpha_polar] = mfradialellipse(initmat, 0, q, 0);

% Plotting
h_fig_plot = figure(3); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', [0.85 0.325 0.098], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Rectangular'); % A warm orange-red
plot(myalpha_radial, falpha_radial, 'Color', [0.494 0.184 0.556], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Radial');   % A deep purple
plot(myalpha_polar, falpha_polar, 'Color', [0.301 0.745 0.933], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Polar');     % A vibrant light blue

hold off; 

xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', fontsz);
legend('Location', 'bestoutside','FontSize',16);
ax = gca;
ax.FontSize = 16; 

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 700, 500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font