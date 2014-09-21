function [atlas, mymaps, p] = align_functional_maps(maps, p, verbose)
%
% function performs groupwise registration among a set of maps
% each map contains
%
% mymaps{i}.Gamma = embedding coordinates
% mymaps{i}.Signals = BOLD signals, tobe compared via correlation
% mymaps{i}.Coords = Coordinates of the points, to be compared via
% Euclidean distance
%
% Transform type:   - orthonormal
%                   - scale_orthonormal
%
% OUTPUT:
% atlas ... structure that holds the registered joint coordinates
% mymaps ... optional output that adds the aligned coordinates to the
% mymaps input structure
%
% Code by Georg Langs and Andy Sweet
% 21.9.2011
%

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

n_subjects = length(mymaps)

if nargin == 1
    p.default = 1;
end
%% Default parameters:
dp.delta = [];
dp.map_dimension = size(mymaps{1}.Gamma,2);
dp.alpha = 0;
%dp.matching_type = 'signal';
%dp.matching_thresh = .6; % for signal
dp.matching_type = 'spatial';
dp.matching_thresh = 40; % for spatial
dp.linear_transform_type = 'orthonormal';
dp.template_subject_index = 1;
dp.groupwise_registration_paradigm = 'all2one';
dp.groupwise_registration_all2all_iterations = 5;
dp.number_of_random_points_in_group_template = ...
    min(size(mymaps{1}.Gamma,1)*3, size(mymaps{1}.Gamma,1)*(length(mymaps)-1));

p = default_parameters(p,dp);
%mymaps.p = p;
%mymaps.Gamma_initial = mymaps.Gamma;
%% Load everything once for all subjects

n_nodes = zeros(n_subjects, 1);
positions = cell(n_subjects, 1);
signals = cell(n_subjects, 1);
funcmaps = cell(n_subjects, 1);
epsilons = cell(n_subjects, 1);
thresholds = cell(n_subjects, 1);

for i=1:n_subjects
    if isfield(mymaps{i},'Signals')
        signals{i} = mymaps{i}.Signals;
    end
    if isfield(mymaps{i},'Coords')
        positions{i} = mymaps{i}.Coords;
    end
    n_nodes(i) = size(mymaps{i}.Gamma,1); % changed from .Coords 5.7.2013 GL
    epsilons{i} = mymaps{i}.Epsilon;
    funcmaps{i} = mymaps{i}.Gamma;
end
clear data
clear map
clear atlas

%% Calculate the linear transforms of the diffusion maps
clear atlas
atlas.Delta = p.delta;


for i=p.template_subject_index
    for j=1:n_subjects
        disp([num2str(i) '-' num2str(j)])
        if ( i ~= j )
            %% 1.) Find correspondences either via spatial proximity of Coords
            % (e.g. if your subjects are registered to MNI space, and you
            % trust the resulting correspondences), or via correlation of
            % BOLD signals (BOLD signals were used for the IPMI2011 results)
            if verbose
                disp('Finding pairwise correspondences ...'); tic;
            end

            if (strcmp(p.matching_type,'spatial')) % use proximity in space as a prior
                [pair_distances, pair_indices_target, pair_indices_source] = sparse_xdistance(positions{i}', positions{j}', p.matching_thresh, verbose);
                epsilon = median(pair_distances);
                pair_weights = exp(-pair_distances/epsilon);
            else % use correlation of signals as a prior
                [pair_corrs, pair_indices_target, pair_indices_source] = sparse_xcorrcoef(signals{i}', signals{j}', p.matching_thresh, verbose);
                pair_weights = exp(pair_corrs / epsilons{i});
            end
            if verbose
                disp(['... done in ', num2str(toc), 's']);
            end

            %% 2.) Calculate the transforms that align map j to the template
            % map via the weighted correspondences
            if verbose
                disp('Linearly aligning correspondences ...'); tic;
            end
            [rotation, translation] = align_points(funcmaps{i}(pair_indices_target,2:end), funcmaps{j}(pair_indices_source, 2:end), pair_weights, p.linear_transform_type);

            % Write the aligned points and the associated cost to the atlas
            % structure
            atlas.LinearRegistrationParams.Rotation{i,j} = eye(p.map_dimension);
            atlas.LinearRegistrationParams.Rotation{i,j}(2:end, 2:end) = rotation(1:p.map_dimension-1,1:p.map_dimension-1);
            atlas.LinearRegistrationParams.Translation{i,j} = zeros(p.map_dimension,1);
            atlas.LinearRegistrationParams.Translation{i,j}(2:end) = translation;
            if verbose
                disp(['... done in ', num2str(toc), 's']);
            end
        end
    end
end

atlas.template_subject_index = p.template_subject_index;
joint_n_nodes = sum(n_nodes);
atlas.Gamma = zeros(joint_n_nodes, p.map_dimension);

%% 3.) Actually apply the transform to the points in the individual maps
atlas.Coords = zeros(size(atlas.Gamma,1),3);
for i=1:n_subjects

    % find indices for this subject in joint
    if ( i == p.template_subject_index )
        atlas.SubjectIndices{i} = 1 : n_nodes(1);
    else
        atlas.SubjectIndices{i} = (sum(n_nodes(1:(i-1))) + 1) : sum(n_nodes(1:i));
    end


       atlas.Coords(atlas.SubjectIndices{i}, :) = positions{i};


    % apply linear transformation to functional map
    % coordinates
    if ( i == p.template_subject_index )
        atlas.Gamma(atlas.SubjectIndices{i}, :) = funcmaps{i};

    else
        atlas.Gamma(atlas.SubjectIndices{i}, :) = ...
            funcmaps{i}*atlas.LinearRegistrationParams.Rotation{p.template_subject_index, i}';
    end
end

for SubjectID = 1:n_subjects
    mymaps{SubjectID}.Gamma_aligned = atlas.Gamma(atlas.SubjectIndices{SubjectID}, :);
end

%% 4.) Perform all to all group-wise registration
if strcmp(p.groupwise_registration_paradigm,'all2all')
    % Currently available for spatial information as basis for
    % registration 22.6.2012

    %% Keeping the all2one gammas:
    atlas.Gamma_all2one = atlas.Gamma;


    for iteration = 1:p.groupwise_registration_all2all_iterations

     %% 4.a)
    random_points = randperm(length(atlas.Gamma));
    random_points = random_points(1:p.number_of_random_points_in_group_template);
    positions_template = atlas.Coords(random_points,:);
    funcmaps_template = atlas.Gamma(random_points,:);

    %% 4.b) recalculate the transformations

        for j=1:n_subjects
            disp([ 'average template - ' num2str(j)])

                %% 1.) Find correspondences either via spatial proximity of Coords
                % (e.g. if your subjects are registered to MNI space, and you
                % trust the resulting correspondences), or via correlation of
                % BOLD signals (BOLD signals were used for the IPMI2011 results)
                if verbose
                    disp('Finding pairwise correspondences ...'); tic;
                end

                if (strcmp(p.matching_type,'spatial')) % use proximity in space as a prior
                    [pair_distances, pair_indices_target, pair_indices_source] = sparse_xdistance(positions_template', positions{j}', p.matching_thresh, verbose);
                    epsilon = median(pair_distances);
                    pair_weights = exp(-pair_distances/epsilon);
                else % use correlation of signals as a prior

                    error('Signal based all2all registration is currently not available') % TODO !!!

                    [pair_corrs, pair_indices_target, pair_indices_source] = sparse_xcorrcoef(signals{i}', signals{j}', p.matching_thresh, verbose);
                    pair_weights = exp(pair_corrs / epsilons{i});
                end
                if verbose
                    disp(['... done in ', num2str(toc), 's']);
                end

                %% 2.) Calculate the transforms that align map j to the template
                % map via the weighted correspondences
                if verbose
                    disp('Linearly aligning correspondences ...'); tic;
                end
                [rotation, translation] = align_points(funcmaps_template(pair_indices_target,2:end), funcmaps{j}(pair_indices_source, 2:end), pair_weights, p.linear_transform_type);

                % Write the aligned points and the associated cost to the atlas
                % structure
                atlas.LinearRegistrationParams.Rotation{1,j} = eye(p.map_dimension);
                atlas.LinearRegistrationParams.Rotation{1,j}(2:end, 2:end) = rotation(1:p.map_dimension-1,1:p.map_dimension-1);
                atlas.LinearRegistrationParams.Translation{1,j} = zeros(p.map_dimension,1);
                atlas.LinearRegistrationParams.Translation{1,j}(2:end) = translation;
                if verbose
                    disp(['... done in ', num2str(toc), 's']);
                end

        end


    %atlas.template_subject_index = p.template_subject_index;
    joint_n_nodes = sum(n_nodes);
    atlas.Gamma_all2all = zeros(joint_n_nodes, p.map_dimension);

    %% 4.c) apply to points

    for i=1:n_subjects

        % find indices for this subject in joint
        if ( i == 1 )
            atlas.SubjectIndices{i} = 1 : n_nodes(1);
        else
            atlas.SubjectIndices{i} = (sum(n_nodes(1:(i-1))) + 1) : sum(n_nodes(1:i));
        end

        % actually apply linear transformation to functional map
        % coordinates

        atlas.Gamma_all2all(atlas.SubjectIndices{i}, :) = ...
            funcmaps{i}*atlas.LinearRegistrationParams.Rotation{p.template_subject_index, i}';

    end

    for SubjectID = 1:n_subjects
        mymaps{SubjectID}.Gamma_all2one_aligned = mymaps{SubjectID}.Gamma_aligned;
        mymaps{SubjectID}.Gamma_aligned = atlas.Gamma(atlas.SubjectIndices{SubjectID}, :);
    end

    atlas.Gamma_all2all_periteration{iteration} = atlas.Gamma_all2all;
    atlas.Gamma = atlas.Gamma_all2all;
    end
end

atlas.Gamma_Info = 'Gamma is equal to Gamma_all2all';
