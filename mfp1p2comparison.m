clear
close all
clc

rng(5); % Set random seed for reproducibility
p1 = 1;
p2 = 0.75;
sz = 2048;
k = 10;
q = -10:0.1:10;

color_rec = [1.0 0.4 0.0]; % Colors
color_rad = [0.0 0.6 0.6]; 
color_pol = [0.7 0.0 1.0];

label_fontsz = 32; % Global label font size
tick_fontsz = 32; % Global tick font size
legend_fontsz = 16; % Global legend font size


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
plot(myalpha_rec, falpha_rec, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular');
plot(myalpha_radial, falpha_radial, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial');
plot(myalpha_polar, falpha_polar, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar'); 

hold off; 

ylim([0 2.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_square_comparison.png', '-dpng', '-r600'); % Save as PNG with 600 DPI

% Plotting Transformed Radial
falpha_radial_transform = falpha_radial + 1;
shift = alpha_theory(101) - myalpha_radial(101);
myalpha_radial_transform = myalpha_radial + shift;

h_fig_plot = figure(2); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular');
plot(myalpha_radial_transform, falpha_radial_transform, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial (Transformed)');
plot(myalpha_polar, falpha_polar, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar');

hold off; 

xlim([1.6 2.5])
ylim([0 2.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_square_comparison_transformed.png', '-dpng', '-r600'); % Save as PNG with 600 DPI
%% Radial MFp1p2 Spectrum Comparison
[initmat, alpha_theory, f_theory] = mfp1p2radial(1024, p1, p2, 8);
[Dq_rec,myalpha_rec,falpha_rec] = mfrectanglebinarized(initmat, 0, q, 0);
[Dq_radial,myalpha_radial,falpha_radial] = mfthetacoordinate(initmat, 0, q, 0);
[Dq_polar,myalpha_polar,falpha_polar] = mfradialellipse(initmat, 0, q, 0);

% Plotting
h_fig_plot = figure(3); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular'); 
plot(myalpha_radial, falpha_radial, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial');   
plot(myalpha_polar, falpha_polar, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar');    

hold off; 

ylim([0 2.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_radial_comparison.png', '-dpng', '-r600'); % Save as PNG with 600 DPI

% Plotting Transformation
falpha_rec_transform = falpha_rec - 1;
falpha_polar_transform = falpha_polar - 1;
shift_polar = myalpha_polar(101) - alpha_theory(101);
shift_rec = myalpha_rec(101) - alpha_theory(101);
myalpha_polar_transform = myalpha_polar - shift;
myalpha_rec_transform = myalpha_rec - shift;

h_fig_plot = figure(4); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec_transform, falpha_rec_transform, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular (Transformed)'); 
plot(myalpha_radial, falpha_radial, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial');  
plot(myalpha_polar_transform, falpha_polar_transform, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar (Transformed)');     

hold off; 

ylim([0 1.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]);
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_radial_comparison_transformed.png', '-dpng', '-r600'); % Save as PNG with 600 DPI
%% Polar MFp1p2 Spectrum Comparison
[initmat, alpha_theory, f_theory] = mfp1p2polar(sz, p1, p2, k);
[Dq_rec,myalpha_rec,falpha_rec] = mfrectanglebinarized(initmat, 0, q, 0);
[Dq_radial,myalpha_radial,falpha_radial] = mfthetacoordinate(initmat, 0, q, 0);
[Dq_polar,myalpha_polar,falpha_polar] = mfradialellipse(initmat, 0, q, 0);

% Plotting
h_fig_plot = figure(5); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular');
plot(myalpha_radial, falpha_radial, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial');
plot(myalpha_polar, falpha_polar, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar');

hold off; 

ylim([0 2.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]); 
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_polar_comparison.png', '-dpng', '-r600'); % Save as PNG with 600 DPI

% Plotting Transformed Radial
falpha_radial_transform = falpha_radial + 1;
shift = alpha_theory(101) - myalpha_radial(101);
myalpha_radial_transform = myalpha_radial + shift;

h_fig_plot = figure(6); 
plot(alpha_theory, f_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on; 

% Using custom RGB colors [R G B]
plot(myalpha_rec, falpha_rec, 'Color', color_rec, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Rectangular');
plot(myalpha_radial_transform, falpha_radial_transform, 'Color', color_rad, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Radial (Transformed)');
plot(myalpha_polar, falpha_polar, 'Color', color_pol, 'LineStyle', ':', 'LineWidth', 1.75, 'DisplayName', 'Polar');

hold off; 

ylim([0 2.1])
legend('Location', 'bestoutside','FontSize',legend_fontsz);
ax = gca;
ax.FontSize = tick_fontsz;
xlabel('$\alpha$', 'Interpreter', 'latex','FontSize', label_fontsz);
ylabel('$f(\alpha)$', 'Interpreter', 'latex','FontSize', label_fontsz);

% Apply general figure settings
set(h_fig_plot, 'Position', [100, 100, 1500, 1500]);
fontname(h_fig_plot, 'Times New Roman'); % Set font

print(h_fig_plot, 'mfp1p2_polar_comparison_transformed.png', '-dpng', '-r600'); % Save as PNG with 600 DPI