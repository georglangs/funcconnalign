classdef SystemException < fcalign.Exception
% a system exception in the fcalign package

methods

function this = SystemException(id, message, varargin)
% creates a system exception in the fcalign package
%
% Arguments
% ---------
% id : string
%   exception identifier
% message : string
%   informative message to present to the user
    parser = inputParser();
    parser.addRequired('id', @(x) validateattributes(x, {'char'}, {'vector'}));
    parser.addRequired('message', @(x) validateattributes(x, {'char'}, {'vector'}));
    parser.parser(id, message);

    systemId = sprintf('system:%s', id);
    this = this@fcalign.Exception(systemId, message, varargin{:});
end

end

end
