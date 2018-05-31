% @file gsTHBSpline.m
%
%    @brief Matlab wrapper for gsTHBSpline class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon, A. Mantzaflaris

classdef gsTHBSpline < handle

    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsTHBSpline(varargin)
            %gsTHBSpline - construct a gsTHBSpline object
            %
            %Usage:
            %  thb = gsTHBSpline( file )
            %
            %Input:
            %  file: char, [1 x numChar].
            %    Name of input file from which to read/construct the
            %    gsTHBSpline.
            %
            %Output:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~(isa(varargin{1},'char')))
                error('Input argument no. 1 should be of type ''char''.')
            elseif (~exist(varargin{1},'file'))
                error('File does not exist: %s.',varargin{1})
            end
            this.objectHandle = mex_gsTHBSpline('constructor', class(varargin{1}), varargin{:});
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsTHBSpline object
            %
            %Usage:
            %  thb.delete()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  (none)
        
            mex_gsTHBSpline('destructor', this.objectHandle);
        end

        % dim - call class method
        function varargout = dim(this, varargin)
            %dim - dimension of the parameter space of a gsTHBSpline object
            %
            %Usage:
            %  val = thb.dim()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the parameter space of the gsTHBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'dim',  varargin{:});
        end
        
        % support - call class method
        function varargout = support(this, varargin)
            %support - support of a gsTHBSpline object
            %
            %Usage:
            %  supp = thb.support()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  supp: double, [1 x 2*d].
            %    Support of the gsTHBSpline, ordered like 
            %      [u1_min, ..., ud_min, u1_max, ..., ud_max]
            %    where d is the parametric dimension of the 
            %    gsTHBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'support',  varargin{:});
        end

        % size - call class method
        function varargout = size(this, varargin)
            %size - size of a gsTHBSpline object
            %
            %Usage:
            %  num = thb.size()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the gsTHBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'size',  varargin{:});
        end

        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsTHBSpline object
            %
            %Usage:
            %  val = thb.eval( pts )
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTHBSpline.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('eval', this.objectHandle, varargin{:});
        end
        
        % deriv - call class method
        function [varargout] = deriv(this, varargin)
            %deriv - evaluate the jacobian of a gsTHBSpline object
            %
            %Usage:
            %  val = thb.deriv( pts )
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTHBSpline.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of the jacobian of all active functions in each of
            %    the specified points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('deriv', this.objectHandle, varargin{:});
        end

        % hess - call class method
        function [varargout] = hess(this, varargin)
            %deriv - evaluate the hessian of a gsTHBSpline object in one
            %direction
            %
            %Usage:
            %  val = thb.hess( pts, dir )
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTHBSpline.
            %  dir: int
            %    Direction of space on which to compute the hessian.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of the hessian matrix in direction dir of all active 
            %    functions in each of the specified points.
            
            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with %d rows.', this.dim())
            end
            if (~isa(varargin{2},'numeric') || ~(mod(varargin{2},1)==0) || varargin{2}>this.dim())
                error('Input argument no. 2 must be an integer smaller than %d.', this.dim())
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('hess', this.objectHandle, varargin{:});
        end
        
        % active - call class method
        function [varargout] = active(this, varargin)
            %active - active functions of a gsTHBSpline object
            %
            %Usage:
            %  act = thb.active( pts )
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the function.
            %
            %Output:
            %  act: double, [numFun x numPts].
            %    Index of active functions in each of the specified points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('active', this.objectHandle, varargin{:});
        end
        
        % basis - call class method
        function [varargout] = basis(this, varargin)
            %basis - returns the gsTHBSplineBasis object linked to a 
            % gsTHBSpline object
            %
            %Usage:
            %  act = thb.basis()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  basis: gsTHBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            basis_ptr = mex_gsTHBSpline('basis', this.objectHandle, varargin{:});
            [varargout{1:nargout}] = gsTHBSplineBasis(basis_ptr);
        end

    end
end
