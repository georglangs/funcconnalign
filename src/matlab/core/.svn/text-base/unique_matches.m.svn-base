function [inds,xs,ys] = unique_matches(x,y,k)

inds = zeros(k,1);
xs = zeros(k,1);
ys = zeros(k,1);

j = 0;
n_matched = 0;

% this seems to be the quickest
nx = numel(x);
duplicate_indices = false(nx,1);
while n_matched <= k
    
    if (~mod(n_matched,250))
        disp(['matched ', num2str(n_matched), ' points']);
    end
    
    if (all(duplicate_indices,1))
        break;
    else
        n_matched = n_matched + 1;
        j = find(~duplicate_indices, 1, 'first');
        inds(n_matched) = j;
        xs(n_matched) = x(j);
        ys(n_matched) = y(j);
        
        duplicate_indices = duplicate_indices | (x == x(j) | y == y(j)); 
    end
end

% % keep looking for next pair, where both values are unique THIS DOESNT WORK
% [ux,uxi] = my_unique(x);
% [uy,uyi] = my_unique(y);
% while n_matched <= k
% 
%     % look for index of next unique pair
%     xj = j + find(uxi(j+1:end), 1, 'first');
%     yj = j + find(uyi(j+1:end), 1, 'first');
%     
%     j = max(xj,yj);
%     
%     % if either are empty, then we're done
%     if(isempty(j)) 
%         break;
%     end
%     
%     % accept if both are unique at this index
%     if ( uxi(j) && uyi(j) )
%         
%         n_matched = n_matched + 1;
%         inds(n_matched) = j;
%         xs(n_matched) = x(j);
%         ys(n_matched) = y(j);
%     end
%     
% end

% while n_matched <= k 
%     
%     % look for next values in x and y that haven't already been added
%     uxis = ~ismember(x(j+1:end), xs(1:n_matched));
%     uyis = ~ismember(y(j+1:end), ys(1:n_matched));
%         
%     ui = (uxis & uyis);
% 
%     % if there is no unique pair remaining, finish
%     if ( isempty(ui) || sum(ui) == 0 )
%         break;
%     else
%     % otherwise add it
%         j = j + find(ui, 1, 'first');
%         
%         n_matched = n_matched + 1;
%         inds(n_matched) = j;
%         xs(n_matched) = x(j);
%         ys(n_matched) = y(j);
%     end
% end

% only keep actual matches
zero_inds = (inds == 0);
inds(zero_inds) = [];
xs(zero_inds) = [];
ys(zero_inds) = [];


%%%%%%%%%%%%%%%%%%%%

% function [ux, ui] = my_unique(x)
% % find unique values and indices in x
% 
% [sx,si] = sort(x, 'ascend');
% 
% ui(si) = ([1; diff(sx)] ~= 0);
% 
% ux = x(ui);
    