% example of functional geometry map building and atlas construction
%
% *** Preliminary release, please don't distribute ***
%
% Details on the method:
%
% G. Langs, D. Lashkari, A. Sweet, Y. Tie, L. Rigolo, A. J. Golby, and
% P. Golland. Learning an atlas of a cognitive process in its functional
% geometry. Inf Process Med Imaging, 22:135?46, 2011.
%
% In this version of the code only orthonormal registration is implemented.
% This will be extended.
%
% Georg Langs
% September 2011
%

clear all
close all

%% Create 2 maps

% The minimal input for map building consists of structures
% data with two fields:
% data.Coords ... spatial coordinates of the points in the brain
% data.BOLD ... BOLD signal at these points
% Typical number of points we used are up to 20000 points.
%

% output messages indicating progress
verbose = true;

% build map for subject 1:
[maps{1}] = create_functional_map(data1)
% add spatial / anatomy coords:
maps{1}.Coords = data1.Coords;

% build map for subject 2:
[maps{2}] = create_functional_map(data2)
% add spatial / anatomy coords:
maps{2}.Coords = data2.Coords;

disp('map building done')

%% Align maps
% orthonormal alignment, init based on proximity in real space
[atlas, mymaps_aligned] = align_group(maps, verbose)


disp('atlas building done')

