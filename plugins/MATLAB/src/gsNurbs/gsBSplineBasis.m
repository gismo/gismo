% @file gsTHBSpline.m
%
%    @brief Matlab wrapper for gsBSplineBasis class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): H.M. Verhelst

classdef gsBSplineBasis < gsBasis

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsBSplineBasis(varargin)
            %gsBSplineBasis - construct a gsBSplineBasis object
            %
            %Usage:
            %  thb = gsBSplineBasis( kv1, kv2 )
            %  thb = gsBSplineBasis( kv1, kv2, kv3 )
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO

            if (~nargin==1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            if (isa(varargin{1},'uint64'))
                this.objectHandle = varargin{1};
            elseif (nargin==1)
                if ( isa(varargin{1},'double') )
                    error('Not implemented');
                    % this.objectHandle = calllib('libgismo','gsBSplineBasis2_create',varargin{:});

                    % this.objectHandle = mex_gsBSplineBasis('constructor', class(varargin{1}), varargin{:});
                elseif ( isa(varargin{1},'gsKnotVector') )
                    this.objectHandle = calllib('libgismo','gsBSplineBasis_create',varargin{1}.ptr());
                else
                    error('Input argument no. 1 should be of type ''double''.')
                end
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Member function implementations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % knots - call class method
        function [varargout] = knots(this, varargin)
            %knots - returns the knot vector of a gsBasis object
            %   of the specified level on the specified direction
            %
            %Usage:
            %  knt = thb.knots( dir )
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %  dir: int, [1 x 1].
            %    index of the direction of space to consider.
            %
            %Output:
            %  knt: double, [1 x numKnots].
            %    Knot vector corresponding to level lev in the direction
            %    dir.

            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}<1 || varargin{1}>this.dim())
                error('Input argument no. 1 must be a strictly positive integer smaller than %d.', this.dim())
            end
            [varargout{1:nargout}] = calllib('libgismo','gsBSplineBasis_knots',this.objectHandle,varargin{1});
        end

        % save - call class method
        function [varargout] = save(this, varargin)
            %save - save a gsTHBSplineBasis object as xml object
            %
            %Usage:
            %  thb.save();
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %   (none - saved into an xml file)

            % if (nargin~=2)
            %     error('Invalid number of input arguments.')
            % end
            % if (~(isa(varargin{1},'char')))
            %     error('Input argument no. 1 should be of type ''char''.')
            % end
            % [varargout{1:nargout}] = mex_gsTHBSplineBasis('save', this.objectHandle, varargin{:});
        end
    end
end
