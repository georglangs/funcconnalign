function [A,B] = generate_testing_data()
% Script to generate two brain fmri signal data for testing functional
% alignment. 
% georg langs
% georg.langs@meduniwien.ac.at
% 22.12.2013
%%


load Hemisphere_Surfaces

addpath ../matlab_packages/matvtk/matlab % <- only for visualization
%%

Surf_LH.values1 = zeros(size(Surf_LH.vertices,1),1);
Surf_LH.values2 = zeros(size(Surf_RH.vertices,1),1);

marker1_LH.R1 = 25>sum((Surf_LH.vertices - repmat([-50 -50 -20],size(Surf_RH.vertices,1),1)).^2,2).^.5
marker1_LH.R2 = 25>sum((Surf_LH.vertices - repmat([-30 -70 60],size(Surf_RH.vertices,1),1)).^2,2).^.5
marker1_LH.R3 = 25>sum((Surf_LH.vertices - repmat([-40 40 40],size(Surf_RH.vertices,1),1)).^2,2).^.5

marker2_LH.R1 = 25>sum((Surf_LH.vertices - repmat([-50 -10 -20],size(Surf_RH.vertices,1),1)).^2,2).^.5
marker2_LH.R2 = 25>sum((Surf_LH.vertices - repmat([-30 -50 50],size(Surf_RH.vertices,1),1)).^2,2).^.5
marker2_LH.R3 = 25>sum((Surf_LH.vertices - repmat([-40 20 30],size(Surf_RH.vertices,1),1)).^2,2).^.5


Surf_LH.values1(find(marker1_LH.R1)) = 1;
Surf_LH.values1(find(marker1_LH.R2)) = 2;
Surf_LH.values1(find(marker1_LH.R3)) = 3;

Surf_LH.values2(find(marker2_LH.R1)) = 1;
Surf_LH.values2(find(marker2_LH.R2)) = 2;
Surf_LH.values2(find(marker2_LH.R3)) = 3;

%
%vtkinit
%vtkplotmesh(Surf_LH.vertices,Surf_LH.faces,Surf_LH.values1)
%vtkplotmesh(Surf_LH.vertices+100,Surf_LH.faces,Surf_LH.values2)


%%



R1sorig = randn(200,1)';
R1s = imfilter(R1sorig,fspecial('gaussian',[1 31],3));
R1s = repmat(R1s,sum(marker1_LH.R1),1);
R1s1 = R1s + .05*randn(size(R1s));

R2sorig = randn(200,1)';
R2s = imfilter(R2sorig,fspecial('gaussian',[1 31],3));
R2s = repmat(R2s,sum(marker1_LH.R2),1);
R2s1 = R2s + .05*randn(size(R2s));

R3sorig = randn(200,1)';
R3s = imfilter(R3sorig,fspecial('gaussian',[1 31],3));
R3s = repmat(R3s,sum(marker1_LH.R3),1);
R3s1 = R3s + .05*randn(size(R3s));


R1s = imfilter(R1sorig,fspecial('gaussian',[1 31],3));
R1s = repmat(R1s,sum(marker2_LH.R1),1);
R1s2 = R1s + .05*randn(size(R1s));

R2s = imfilter(R2sorig,fspecial('gaussian',[1 31],3));
R2s = repmat(R2s,sum(marker2_LH.R2),1);
R2s2 = R2s + .05*randn(size(R2s));

R3s = imfilter(R3sorig,fspecial('gaussian',[1 31],3));
R3s = repmat(R3s,sum(marker2_LH.R3),1);
R3s2 = R3s + .05*randn(size(R3s));

hold off
plot(R1s1','r')
hold on
plot(R2s1','g')
plot(R3s1','b')

%%
BOLD1 = .3*randn(size(Surf_LH.vertices,1),200);
BOLD1(find(marker1_LH.R1),:) = R1s1;
BOLD1(find(marker1_LH.R2),:) = R2s1;
BOLD1(find(marker1_LH.R3),:) = R3s1;

BOLD2 = .3*randn(size(Surf_LH.vertices,1),200);
BOLD2(find(marker2_LH.R1),:) = R1s2;
BOLD2(find(marker2_LH.R2),:) = R2s2;
BOLD2(find(marker2_LH.R3),:) = R3s2;



A.BOLD = BOLD1;%rand(400,100);
B.BOLD = BOLD2;%rand(400,100);
A.Coords = Surf_LH.vertices;
B.Coords = Surf_LH.vertices;
A.affinity_matrix = corrcoef(A.BOLD');
B.affinity_matrix = corrcoef(B.BOLD');
A.marker = marker1_LH;
B.marker = marker2_LH;




