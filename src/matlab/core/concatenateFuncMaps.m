function [jointmap] = concatenateFuncMaps(map)
%
% function to concatenate the coordinats of multiple functional maps
%
% input: 
% map{1}.Phi, map{2}.Phi ....
%
% output:
% jointmap.Phi
% jointmap.Pointer{1} etc
%
% langs@csail.mit.edu
% 1.12.2010
%

jointmap.Phi = [];

for count = 1:length(map)
    jointmap.Pointer{count} = [size(jointmap.Phi,1)+1:size(jointmap.Phi,1)+size(map{count}.Phi,1)];
    jointmap.Phi = [jointmap.Phi;map{count}.Phi];
end
    