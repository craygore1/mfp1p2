clear
close all
clc

rng(5); % Set random seed for reproducibility
%% Square MFp1p2

p1 = 1;
p2 = 0.75;
sz = 2048;

% Values of k to plot
k_values = [1, 2, 6, 10];
num_k_values = length(k_values);

for idx = 1:num_k_values
    k_val = k_values(idx);
    
    initmat = mfp1p2gen(sz, p1, p2, k_val);
    
    h_fig = figure; % Get handle to the current figure
    imshow(initmat);

    filename = sprintf('Square_MFp1p2_k%d.png', k_val);
    
    print(h_fig, filename, '-dpng', '-r600'); % Save as PNG with 600 DPI
    
    close(h_fig);
end

%% Radial MFp1p2
p1 = 1;
p2 = 0.75;
sz = 2048;

% Values of k to plot
k_values = [1, 2, 4, 8];
num_k_values = length(k_values);

for idx = 1:num_k_values
    k_val = k_values(idx);
    
    initmat = mfp1p2radial(sz, p1, p2, k_val);
    
    h_fig = figure; % Get handle to the current figure
    imshow(initmat);
    
    filename = sprintf('Radial_MFp1p2_k%d.png', k_val);
    
    print(h_fig, filename, '-dpng', '-r600'); % Save as PNG with 600 DPI
    
    close(h_fig);
end

%% Polar MFp1p2
p1 = 1;
p2 = 0.75;
sz = 2048;

% Values of k to plot
k_values = [2, 4, 6, 10];
num_k_values = length(k_values);

for idx = 1:num_k_values
    k_val = k_values(idx);
    
    initmat = mfp1p2polar(sz, p1, p2, k_val);
    
    h_fig = figure; % Get handle to the current figure
    imshow(initmat);
    
    filename = sprintf('Polar_MFp1p2_k%d.png', k_val);
    
    print(h_fig, filename, '-dpng', '-r600'); % Save as PNG with 600 DPI
    
    close(h_fig);
end