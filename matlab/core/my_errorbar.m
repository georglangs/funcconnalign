function my_errorbar(x, means, mean_line_width, stds, std_line_width, color)

% plot the means
plot(x, means, 'LineWidth', mean_line_width, 'color', color);

% plot vertical line for std of each mean
hold on;
for i=1:numel(x)
    X = [x(i), x(i)];
    Y = [means(i)-stds(i), means(i) + stds(i)];
    line(X, Y, 'LineWidth', std_line_width, 'color', color);
end
