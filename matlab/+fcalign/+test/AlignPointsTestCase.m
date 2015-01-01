classdef AlignPointsTestCase < fcalign.test.TestCase
% tests different usages of aligning multiple functional maps

methods (Test)

function test_default(this)
% tests default usage of aligning multiple functional maps

    % define parameters
    n = 16;
    p = 3;

    % create some data points
    points = this.create_data(n, p);

    % define a transform and transform the points
    expected_transform = fcalign.create_transform(p);
    transformed_points = fcalign.transform_points(points, expected_transform);

    % estimate the transform
    actual_transform = fcalign.align_points(points, transformed_points);

    % verify that the transform
    this.verify_equal_transform(actual_transform, expected);
end

end

end
