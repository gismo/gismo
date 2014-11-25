% MATLAB gsTHBSplineBasis class wrapper for the C++ gsTHBSplineBasis class

% Author: Peter Noertoft.

classdef gsTHBSplineBasis < handle

    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

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
            if (~(isa(varargin{1},'char')))
                error('Input argument no. 1 should be of type ''char''.')
            elseif (~exist(varargin{1},'file'))
                error('File does not exist: %s.',varargin{1})
            end
            this.objectHandle = mex_gsTHBSplineBasis('constructor', class(varargin{1}), varargin{:});
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

        % dim - call class method
        function varargout = dim(this, varargin)
            %dim - dimension of the parameter space of a gsTHBSplineBasis object
            %
            %Usage:
            %  val = thb.dim()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the parameter space of the gsTHBSplineBasis.
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'dim',  varargin{:});
        end
        
        % numElements - call class method
        function varargout = numElements(this, varargin)
            %numElements - number of elements of a gsTHBSplineBasis object
            %
            %Usage:
            %  num = thb.numElements()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Number of elements of the gsTHBSplineBasis.
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'numElements',  varargin{:});
        end
        
        % support - call class method
        function varargout = support(this, varargin)
            %support - support of a gsTHBSplineBasis object
            %
            %Usage:
            %  supp = thb.support()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  supp: double, [1 x 2*d].
            %    Support of the gsTHBSplineBasis, ordered like 
            %      [u1_min, ..., ud_min, u1_max, ..., ud_max]
            %    where d is the parametric dimennsion of the 
            %    gsTHBSplineBasis.

            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'support',  varargin{:});
        end

        % size - call class method
        function varargout = size(this, varargin)
            %size - size of a gsTHBSplineBasis object
            %
            %Usage:
            %  num = thb.size()
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the gsTHBSplineBasis.
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'size',  varargin{:});
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
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'treeSize',  varargin{:});
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
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('accessor', this.objectHandle, 'treeLeafSize',  varargin{:});
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
        
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('treePrintLeaves', this.objectHandle, varargin{:});
        end

        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsTHBSplineBasis object
            %
            %Usage:
            %  val = thb.eval( pts )
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTHBSplineBasis.
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
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('eval', this.objectHandle, varargin{:});
        end

        % evalSingle - call class method
        function [varargout] = evalSingle(this, varargin)
            %evalSingle - evaluate a single function in a gsTHBSplineBasis object
            %
            %Usage:
            %  valSingle = thb.evalSingle( fun, pts )
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
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
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}))
                error('Input argument no. 1 must be a numeric scalar.')
            elseif (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
                error('Input argument no. 2 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('evalSingle', this.objectHandle, varargin{:});
        end

        % active - call class method
        function [varargout] = active(this, varargin)
            %active - active functions of a gsTHBSplineBasis object
            %
            %Usage:
            %  act = thb.active( pts )
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
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
            [varargout{1:nargout}] = mex_gsTHBSplineBasis('active', this.objectHandle, varargin{:});
        end

    end
end