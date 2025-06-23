clear
close all
clc

rng(5); % Set random seed for reproducibility
p1 = 1;
p2 = 0.75;
sz = 2048;

% Values of k to plot
k_values = [1, 2, 6, 10];

num_k_values = length(k_values);

num_rows = ceil(sqrt(num_k_values));
num_cols = ceil(num_k_values / num_rows);

t = tiledlayout(num_rows, num_cols, ...
    'TileSpacing', 'compact', 'Padding', 'tight');     

for idx = 1:num_k_values
    k_val = k_values(idx);
    
    initmat = mfp1p2gen(sz, p1, p2, k_val);

    nexttile; % Moves to the next tile in the layout
    imshow(initmat);
    xlabel(['k = ' num2str(k_val)], 'FontSize', 24);
    
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
end

set(gcf, 'Position', [100, 100, 1000, 1000])