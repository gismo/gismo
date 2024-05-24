% @file gsKnotVector.m
%
%    @brief Matlab wrapper for gsKnotVector class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): H.M. Verhelst

classdef gsKnotVector < handle
% try: Custom Display Interface ? override disp ?
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsKnotVector(varargin)
            %gsKnotVector - construct a gsKnotVector object
            %
            %Usage:
            %  gsKv = gsKnotVector( deg,kv )
            %
            %Input:
            %  deg: int, [1x1]
            %    Degree of the knot vector
            %  kv: real, [numKnots x 1].
            %    Row vector representing the knots
            %
            %Output:
            %  gsKv: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            % HV: WHAT DOES THIS ONE DO??
            if (isa(varargin{1},'uint64'))
                this.objectHandle = varargin{1};
            else
                if (~(isa(varargin{1},'double')))
                    error('Input argument no. 2 should be of type ''double''.')
                end
                this.objectHandle = calllib('libgismo','gsKnotVector_create',varargin{1},length(varargin{1}));
            end
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsKnotVector object
            %
            %Usage:
            %  thb.delete()
            %
            %Input:
            %  thb: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            %
            %Output:
            %  (none)

            calllib('libgismo','gsKnotVector_delete',this.ptr);
        end

         % This function returns the address of the C++ pointer
         function varargout = ptr(this, varargin)
             if (nargin~=1 || nargout>1)
                 error('Invalid number of input and/or output arguments.')
             end
             varargout{1} = this.objectHandle;
         end

         function disp(this, var_name)
            calllib('libgismo','gsKnotVector_print',this.ptr);
             %a = mex_gsKnotVector('get',this.objectHandle,[]);
             %fmt=['gsknotVector {' repmat(' %1.0f',1,numel(a)) ' }\n'];
             %fprintf(fmt,a);
        end

        % degree - call class method
        function [varargout] = degree(this, varargin)
            %degree - the degree for a specified direction of a 
            %   gsKnotVector object
            %
            %Usage:
            %  deg = thb.degree( dir )
            %
            %Input:
            %  thb: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            %  dir: int, [1 x 1].
            %    Direction of space for which we want to know the degree of 
            %    the gsKnotVector.
            %
            %Output:
            %  deg: double, [1 x 1].
            %    Degree of the gsKnotVector object in direction dir.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            % TODO
            % [varargout{1:nargout}] = calllib('libgismo','gsKnotVector_degree',this.ptr);
        end
    end
end
