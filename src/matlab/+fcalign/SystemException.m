classdef SystemException < fcalign.Exception
    % A system exception in the fcalign package.

methods

function this = SystemException(id, message, varargin)
    % Creates a system exception in the fcalign package.
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

    systemId = sprintf('system:%s', id);
    this = this@fcalign.Exception(systemId, message, varargin{:});
end

end

end
