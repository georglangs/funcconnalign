classdef Exception < MException
    % An exception in the fcalign package.

methods

function this = Exception(id, message, varargin)
    % Creates an exception in the fcalign package.
    %
    % Arguments
    % ---------
    % id : str
    %   Exception identifier.
    % message : str
    %   Informative message to present to the user.

    parser = inputParser();
    parser.addRequired('id', @(x) validateattributes(x, {'char'}, {'vector'}));
    parser.addRequired('message', @(x) validateattributes(x, {'char'}, {'vector'}));
    parser.parse(id, message);

    fullId = sprintf('fcalign:%s', id);
    this = this@MException(fullId, message, varargin{:});
end

end

end
