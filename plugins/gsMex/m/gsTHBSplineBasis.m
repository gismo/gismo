% @file gsTHBSpline.m
%
%    @brief Matlab wrapper for gsTHBSplineBasis class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon, P. Noertoft

classdef gsTHBSplineBasis < gsBasis

    % properties (SetAccess = private, Hidden = true)
    %     objectHandle; % Handle to the underlying C++ class instance
    % end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsTHBSplineBasis(varargin)
            %gsTHBSplineBasis - construct a gsTHBSplineBasis object
            %
            %Usage:
            %  thb = gsTHBSplineBasis( file )
            %
            %Input:
            %  file: char, [1 x numChar].
            %    Name of input file from which to read/construct the
            %    gsTHBSplineBasis.
            %
            %Output:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            if (isa(varargin{1},'uint64'))
                this.objectHandle = varargin{1};
            else
                if ( isa(varargin{1},'char') ) % construct using XML filename string
                    % if exist(varargin{1},'file') then read from file..
                    this.objectHandle = mex_gsTHBSplineBasis('constructor', class(varargin{1}), varargin{:});
                elseif ( isa(varargin{1},'gsTensorBSplineBasis') )
                    this.objectHandle = mex_gsTHBSplineBasis('constructor', class(varargin{1}), varargin{1}.ptr());
                else
                    error('Input argument no. 1 should be of type char, file or gsTensorBSplineBasis.')
                end
            end
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsTHBSplineBasis object
            %
            %Usage:
            %  thb.delete()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  (none)
        
            mex_gsTHBSplineBasis('destructor', this.objectHandle);
        end

         % This function returns the address of the C++ pointer
         function varargout = ptr(this, varargin)
             if (nargin~=1 || nargout>1)
                 error('Invalid number of input and/or output arguments.')
             end
             varargout{1} = this.objectHandle;
         end

        % treeSize - call class method
        function varargout = treeSize(this, varargin)
            %treeSize - size of the tree of a gsTHBSplineBasis object
            %
            %Usage:
            %  num = thb.treeSize()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the tree of the gsTHBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('treeSize', this.objectHandle, varargin{:});
        end
        
        % treeLeafSize - call class method
        function varargout = treeLeafSize(this, varargin)
            %treeLeafSize - size of the leaf in the tree of a gsTHBSplineBasis object
            %
            %Usage:
            %  num = thb.treeLeafSize()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the leaf in the tree of the gsTHBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('treeLeafSize', this.objectHandle,  varargin{:});
        end

        % maxLevel - call class method
        function varargout = maxLevel(this, varargin)
            %maxLevel - maximum level of a gsTHBSplineBasis object
            %
            %Usage:
            %  lev = thb.maxLevel()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  lev: double, [1 x 1].
            %    Maximum level present in the hierarchy of the
            %    gsTHBSplineBasis object. 
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('maxLevel', this.objectHandle,  varargin{:});
        end
        
        % treePrintLeaves - call class method
        function varargout = treePrintLeaves(this, varargin)
            %treePrintLeaves - print the leaves in the tree of a gsTHBSplineBasis object
            %
            %Usage:
            %  thb.treePrintLeaves()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  (none - outputs to the screen).
            
            if (nargin~=1 || nargout>0)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('treePrintLeaves', this.objectHandle, varargin{:});
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
            
            if (nargin~=2)
                error('Invalid number of input arguments.')
            end
            if (~(isa(varargin{1},'char')))
                error('Input argument no. 1 should be of type ''char''.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('save', this.objectHandle, varargin{:});
        end
        
        % knots - call class method
        function [varargout] = knots(this, varargin)
            %knots - returns the knot vector of a gsTHBSplineBasis object
            %   of the specified level on the specified direction
            %
            %Usage:
            %  knt = thb.knots( lev, dir )
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %  lev: int, [1 x 1].
            %    Index of the level of the hierarchy to consider.
            %  dir: int, [1 x 1].
            %    index of the direction of space to consider. 
            %
            %Output:
            %  knt: double, [1 x numKnots].
            %    Knot vector corresponding to level lev in the direction
            %    dir. 
            
            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ~(mod(varargin{1},1)==0) || varargin{1}<1)
                error('Input argument no. 1 must be a strictly positive integer.')
            elseif (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || ...
                    ~(mod(varargin{2},1)==0) || varargin{2}<1 || varargin{2}>this.dim())
                error('Input argument no. 2 must be a strictly positive integer smaller than %d.', this.dim())
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('knots', this.objectHandle, varargin{:});
        end
    end
end