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
%    Author(s): H.M.Verhelst

classdef gsBasis < gsFunctionSet

    properties (SetAccess = public, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

        % dim - call class method
        function varargout = dim(this, varargin)
            %dim - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsBasis_dim',this.objectHandle);
        end
        
        % numElements - call class method
        function varargout = numElements(this, varargin)
            %numElements - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsBasis_numElements',this.objectHandle);
        end

        % uniformRefine - call class method
        function uniformRefine(this, varargin)
            %uniformRefine - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO

            if (~(nargin>=1 && nargin <= 4) || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            if (nargin<4) dir = -1; else dir = varargin{3}; end
            if (nargin<3) mul =  1; else mul = varargin{2}; end
            if (nargin<2) num =  1; else num = varargin{1}; end

            calllib('libgismo','gsBasis_uniformRefine',this.objectHandle,num,mul,dir);
        end

        % refineElements - call class method
        function refineElements(this, varargin)
            %refineElements - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO

            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end

            boxes = int32(varargin{1});
            len = length(boxes);
            if (~isa(boxes,'integer') || ~isvector(boxes) || ~isequal(mod(len,2*this.domainDim()+1),0))
                error('Input argument no. 1 must be an integer vector and with 2*d+1 rows or cols.')
            end
            calllib('libgismo','gsBasis_refineElements',this.objectHandle,varargin{1},len);
        end
        
        % refineElements - call class method
        function refine(this, varargin)
            %refineElements - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO

            if (~(nargin==2 || nargin==3) || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end

            if (nargin<3) refExt = 0; else dir = varargin{2}; end

            len = length(varargin{1});
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}))
                error('Input argument no. 1 must be a numeric matrix.')
            end

            assert(isequal(size(varargin{1},1),this.domainDim()) && isequal(size(varargin{1},2),2),"Boxes must be domainDim() by 2")

            boxes = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            calllib('libgismo','gsBasis_refine',this.objectHandle,boxes.ptr(),refExt);
        end

        % support - call class method
        function varargout = support(this, varargin)
            %support - TODO
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO
            
            if (~(nargin==1 || nargin==2) || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (nargin==1)
                [varargout{1:nargout}] = calllib('libgismo','gsFunctionSet_support',this.objectHandle);
            elseif (nargin==2)
                [varargout{1:nargout}] = calllib('libgismo','gsBasis_support',this.objectHandle, varargin{:});
            end
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
            [varargout{1:nargout}] = calllib('libgismo','gsBasis_size',this.objectHandle);
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
            [varargout{1:nargout}] = calllib('libgismo','gsBasis_degree',this.objectHandle, varargin{:});
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

            assert(size(varargin{1},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrixInt();
            calllib('libgismo','gsBasis_active_into',this.objectHandle,uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
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

            if (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{2},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{2},1),size(varargin{2},2),varargin{2});
            result = EigenMatrix();
            calllib('libgismo','gsBasis_evalSingle_into',this.objectHandle,varargin{1},uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

        % derivSingle - call class method
        function [varargout] = derivSingle(this, varargin)
            %derivSingle - evaluate a single function in a gsBasis object
            %
            %Usage:
            %  valSingle = thb.derivSingle( fun, pts )
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

            if (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{2},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{2},1),size(varargin{2},2),varargin{2});
            result = EigenMatrix();
            calllib('libgismo','gsBasis_derivSingle_into',this.objectHandle,varargin{1},uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

        % deriv2Single - call class method
        function [varargout] = deriv2Single(this, varargin)
            %deriv2Single - evaluate a single function in a gsBasis object
            %
            %Usage:
            %  valSingle = thb.deriv2Single( fun, pts )
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

            if (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{2},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{2},1),size(varargin{2},2),varargin{2});
            result = EigenMatrix();
            calllib('libgismo','gsBasis_deriv2Single_into',this.objectHandle,varargin{1},uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

    end
end
