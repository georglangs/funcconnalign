function [y Flag] = logbz(x)

% [y Flag] = logbz(x)
%
% This is a logarithm function except that it is bounded below by -200 when
% x < 1e-200 in which case Flag = 1
% 

Flag = (x < 1e-200);
y = log( ~Flag .* x + Flag .* 1e-200 ); 