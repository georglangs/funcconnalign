function trans_pts = fsl_trans_image_points(dest, src, trans, pts, trans_type)
% transforms points using FSL. This assumes all points are specific in real world mm coordinates. It writes temporary text files in the current file and deletes them before finishing.
%
% INPUTS
% dest          filename of destination image
% src           filename of source image
% trans         filename of transformat
% pts           set of pts (Nxd : X d-dim points) should d=3 only?
% trans_type    'linear': linear transformation matrix in format output by FLIRT (default)
%               'warp': warp field in format output bY FNIRT
%
% OUTPUTS
% trans_pts     set of transformed points
% 

fsl_dir='/csail/vision-polina4/shared_software/fsl';

if (nargin < 5)
    trans_type = 'linear';
end
MyToken = round(rand*1000000)
pts_file = ['pts_src_temp' num2str(MyToken) '.txt'];
save(pts_file, 'pts', '-ascii');


trans_pts_file = ['pts_src_trans_temp_' num2str(MyToken) '.txt'];
trans_pts_file2 = ['pts_src_trans_temp2_' num2str(MyToken) '.txt'];

if (strcmp(trans_type, 'warp'))
    unix([fsl_dir, '/bin/img2imgcoord -mm -src ', src, ' -dest ', dest, ' -warp ', trans, ' ', pts_file, ' > ', trans_pts_file]);
else
    unix([fsl_dir, '/bin/img2imgcoord -mm -src ', src, ' -dest ', dest, ' -xfm ', trans, ' ', pts_file, ' > ', trans_pts_file]);
end

% get rid of first line
unix(['more +2 ', trans_pts_file, ' > ', trans_pts_file2]); 

% load from text file and return points
trans_pts = load(trans_pts_file2);

% get rid of any repeats at ends (why does this happen?)
n_pts = size(pts,1);
n_trans_pts = size(trans_pts,1);
if (n_pts < n_trans_pts)
    trans_pts = trans_pts(1:n_pts,:);
end

% remove temporary files
unix(['rm ', pts_file]);
unix(['rm ', trans_pts_file]);
unix(['rm ', trans_pts_file2]);
