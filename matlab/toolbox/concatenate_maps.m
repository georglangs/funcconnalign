function map0 = concatenate_maps(map1,map2)
% function to concatenate functional connectivity maps
% georg langs
% langs@csail.mit.edu
% 26.7.2012
%

map0.NodeIndices = [map1.NodeIndices;map2.NodeIndices];

map0.concatenation_split = [size(map1.Gamma,1)];

if map1.Epsilon == map2.Epsilon
    map0.Epsilon = map1.Epsilon;
else
    map0.Epsilon = [map1.Epsilon map2.Epsilon];
end

if map1.CorrelationThreshold == map2.CorrelationThreshold
    map0.CorrelationThreshold = map1.CorrelationThreshold;
else
    map0.CorrelationThreshold = [map1.CorrelationThreshold map2.CorrelationThreshold];
end

if strcmp(map1.ThresholdType,map2.ThresholdType)
    map0.ThresholdType = map1.ThresholdType;
else
    error('different threshold types, not concatenating')
    %map0.CorrelationThreshold = [map1.CorrelationThreshold map2.CorrelationThreshold];
end

if map1.MinimumNodeDegree == map2.MinimumNodeDegree
    map0.MinimumNodeDegree = map1.MinimumNodeDegree;
else
    error('different MinimumNodeDegree, not concatenating')
    %map0.CorrelationThreshold = [map1.CorrelationThreshold map2.CorrelationThreshold];
end

if strcmp(map1.SamplingType,map2.SamplingType)
    map0.SamplingType = map1.SamplingType;
else
    error('different SamplingType, not concatenating')
    %map0.CorrelationThreshold = [map1.CorrelationThreshold map2.CorrelationThreshold];
end

if map1.Alpha == map2.Alpha
    map0.Alpha = map1.Alpha;
else
    error('different Alpha, not concatenating')
    %map0.CorrelationThreshold = [map1.CorrelationThreshold map2.CorrelationThreshold];
end

map0.GraphDegrees = [map1.GraphDegrees;map2.GraphDegrees];
map0.Degrees = [map1.Degrees;map2.Degrees];
map0.MsEigenvalues = [map1.MsEigenvalues;map2.MsEigenvalues];
map0.Gamma = [map1.Gamma;map2.Gamma];


if map1.DiffusionTime == map2.DiffusionTime
    map0.DiffusionTime = map1.DiffusionTime;
else
    map0.DiffusionTime = [map1.DiffusionTime map2.DiffusionTime];
end

if map1.DeterminedDelta == map2.DeterminedDelta
    map0.DeterminedDelta = map1.DeterminedDelta;
else
    map0.DeterminedDelta = [map1.DeterminedDelta map2.DeterminedDelta];
end

if map1.Delta == map2.Delta
    map0.Delta = map1.Delta;
else
    map0.Delta = [map1.Delta map2.Delta];
end


map0.MapDimension = map1.MapDimension;
map0.Coords = [map1.Coords;map2.Coords];
if isfield(map1,'Signals') & isfield(map2,'Signals')
map0.Signals = [map1.Signals;map2.Signals];
end

