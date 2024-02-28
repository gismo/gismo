% @file gsBasis.m
%
%    @brief Matlab wrapper for gsBasis class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon, P. Noertoft

classdef gsBasis < handle

    properties (SetAccess = public, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

         % This function returns the address of the C++ pointer
         function varargout = ptr(this, varargin)
             if (nargin~=1 || nargout>1)
                 error('Invalid number of input and/or output arguments.')
             end
             varargout{1} = this.objectHandle;
         end

        % dim - call class method
        function varargout = dim(this, varargin)
            %dim - dimension of the parameter space of a gsBasis object
            %
            %Usage:
            %  val = thb.dim()
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the parameter space of the gsBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsBasis('dim', this.objectHandle,  varargin{:});
            %mex_gsBasis('dim', this.objectHandle, varargin{:});
        end
        
        % numElements - call class method
        function varargout = numElements(this, varargin)
            %numElements - number of elements of a gsBasis object
            %
            %Usage:
            %  num = thb.numElements()
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Number of elements of the gsBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsBasis('numElements', this.objectHandle, varargin{:});
        end
        
        % support - call class method
        function varargout = support(this, varargin)
            %support - support of a gsBasis object
            %
            %Usage:
            %  supp = thb.support()
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %
            %Output:
            %  supp: double, [1 x 2*d].
            %    Support of the gsBasis, ordered like
            %      [u1_min, ..., ud_min, u1_max, ..., ud_max]
            %    where d is the parametric dimennsion of the 
            %    gsBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsBasis('support', this.objectHandle,  varargin{:});
        end

        % size - call class method
        function varargout = size(this, varargin)
            %size - size of a gsBasis object
            %
            %Usage:
            %  num = thb.size()
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the gsBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsBasis('size', this.objectHandle,  varargin{:});
        end

        % degree - call class method
        function [varargout] = degree(this, varargin)
            %degree - the degree for a specified direction of a 
            %   gsBasis object
            %
            %Usage:
            %  deg = thb.degree( dir )
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %  dir: int, [1 x 1].
            %    Direction of space for which we want to know the degree of 
            %    the gsBasis.
            %
            %Output:
            %  deg: double, [1 x 1].
            %    Degree of the gsBasis object in direction dir.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}>this.dim())
                error('Input argument must be an integer less than %d.', this.dim())
            end
            [varargout{1:nargout}] = mex_gsBasis('degree', this.objectHandle, varargin{:});
        end
        
        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsBasis object
            %
            %Usage:
            %  val = thb.eval( pts )
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsBasis.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.
            disp(this.dim());
            disp(size(varargin{1},1));
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsBasis('eval', this.objectHandle, varargin{:});
        end

        % evalSingle - call class method
        function [varargout] = evalSingle(this, varargin)
            %evalSingle - evaluate a single function in a gsBasis object
            %
            %Usage:
            %  valSingle = thb.evalSingle( fun, pts )
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
            %  fun: double, [1 x 1].
            %    Index of function to evaluate.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the function.
            %
            %Output:
            %  val: double, [1 x numPts].
            %    Value of the specified function in each of the specified
            %    points.
            
            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}<1)
                error('Input argument no. 1 must be an strictly positive integer.')
            elseif (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
                error('Input argument no. 2 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsBasis('evalSingle', this.objectHandle, varargin{:});
        end
        % active - call class method
        function [varargout] = active(this, varargin)
            %active - active functions of a gsBasis object
            %
            %Usage:
            %  act = thb.active( pts )
            %
            %Input:
            %  thb: gsBasis, [1 x 1].
            %    The gsBasis object.
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
            [varargout{1:nargout}] = mex_gsBasis('active', this.objectHandle, varargin{:});
        end

    end
end