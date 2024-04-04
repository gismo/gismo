classdef Geometry < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = public, Hidden = true)
        ptr
    end

    methods
        function this = Geometry(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~(isa(varargin{1},'string')) && ~(isa(varargin{1},'char')))
                error('Input argument no. 1 should be of type ''char''.')
            elseif (~exist(varargin{1},'file'))
                error('File does not exist: %s.',varargin{1})
            end
            if (isa(varargin{1},'string'))
                name = convertStringsToChars(varargin{1});
                this.ptr = calllib('libgismo','gsCReadFile',name);
            elseif (isa(varargin{1},'char'))
                this.ptr = calllib('libgismo','gsCReadFile',varargin{1});
            end

        end

        function delete(this)
            calllib('libgismo','gsFunctionSet_delete',this.ptr);
        end
    end
end