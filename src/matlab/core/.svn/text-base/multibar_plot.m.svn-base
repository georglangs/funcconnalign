function multibar_plot(x, ys, line_width, color, opacities)
% bar plot, where the number of bars at each point on the x axis can vary

n =  numel(x);

% find the numbers of the bars at each point
ks = zeros(n,1);
for i=1:n
    ks(i) = numel(ys{i});
end

% determine the width allocated for each set of bars based on the maximum number of bars we need to plot at a single point
max_k = max(ks);
if length(x)>1
x_width = x(2)-x(1);
else
    x_width = 1;
end
bar_sep_dist = (0.75*x_width)/max_k;

% colormap_lut = sqrt(colormap(hot(100)));
colormap_lut = colormap(jet(100));

% axes, hold on;
for i=1:n
   if not(isempty(ys{i})) 
    % determine the xs (i.e. where to plot the bars)
    i_width = ks(i)*bar_sep_dist;
    xs = linspace(x(i) - i_width/2, x(i) + i_width/2, ks(i));
    
    for j=1:ks(i)
        X = [xs(j), xs(j)];
        Y = [0, ys{i}(j)];
        if (nargin < 5)
            line(X, Y, 'LineWidth', line_width, 'color', color);
        else
            col_index = floor(opacities{i}(j)*100);
            col_index = max(col_index,1);
            line(X, Y, 'LineWidth', line_width, 'color', colormap_lut(col_index,:));
        end
    end
   end
end
