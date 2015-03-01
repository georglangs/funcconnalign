classdef AlignPointsTestCase < fcalign.test.TestCase
    % Tests different usages of aligning point sets.

methods (Test)

function test_default(this)
    % Tests default usage of aligning two point sets.

    % define parameters
    n = 3;
    s = 16;

    % create some data points
    points = this.create_data(n, s);

    % define a transform and transform the points
    expected_transform = this.create_transform(n);
    transformed_points = fcalign.transform_points(points, expected_transform);

    % estimate the transform
    actual_transform = fcalign.align_points(points', transformed_points');

    % verify that the transform
    this.verify_equal_transform(actual_transform, expected_transform);
end

end

end
