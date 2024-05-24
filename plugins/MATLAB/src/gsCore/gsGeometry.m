% @file gsGeometry.m
%
%    @brief Matlab wrapper for gsGeometry class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): H.M.Verhelst

classdef gsGeometry < gsFunctionSet

    properties (SetAccess = public, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

        % basis - call class method
        function varargout = basis(this, varargin)
            %basis - returns a pointer to the basis
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
            [varargout{1:nargout}] = calllib('libgismo','gsGeometry_basis',this.objectHandle);
        end

        % coefs - call class method
        function varargout = coefs(this, varargin)
            %coefs - returns the coefs
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
            result = EigenMatrix();
            calllib('libgismo','gsGeometry_coefs_into',this.objectHandle,result);
            [varargout{1:nargout}] = result.asMatrix();
        end

        % uniformRefine - call class method
        function varargout = uniformRefine(this, varargin)
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

            if (~(nargin>=0 && nargin <= 3) || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (nargin<4) dir = -1; else dir = varargin{3}; end
            if (nargin<3) mul =  1; else mul = varargin{2}; end
            if (nargin<2) num =  1; else num = varargin{1}; end

            calllib('libgismo','gsGeometry_uniformRefine',this.objectHandle,num,mul,dir);
        end

        % refineElements - call class method
        function varargout = refineElements(this, varargin)
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

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end

            len = length(varargin{1});
            if (~isa(varargin{1},'numeric') || ~isvector(varargin{1}) || ~isequal(mod(len,2*this.domainDim()+1),0))
                error('Input argument no. 1 must be a numeric vector and with 2*d+1 rows or cols.')
            end
            boxes = EigenVector(len,varargin{1});
            calllib('libgismo','gsGeometry_refineElements',this.objectHandle,boxes,len);
        end

        % evalSingle - call class method
        function [varargout] = normal(this, varargin)
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

            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end

            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.domainDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{1},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrix();
            calllib('libgismo','gsGeometry_normal_into',this.objectHandle,uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

        % evalSingle - call class method
        function [varargout] = closest(this, varargin)
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

            if (~(nargin==2 || nargin==3) || nargout>2)
                error('Invalid number of input and/or output arguments.')
            end

            if (~isa(varargin{1},'numeric') || ~isvector(varargin{1}) || ~isequal(length(varargin{1}),this.targetDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with targetDim() rows.')
            end

            if (nargin<3) accuracy = 1e-6; else accuracy = varargin{2}; end

            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrix();
            varargout{1} = calllib('libgismo','gsGeometry_closestPointTo',this.objectHandle,uu.ptr(),result.ptr(),accuracy);
            varargout{2} = result.asMatrix();
        end

        % evalSingle - call class method
        function [varargout] = invertPoints(this, varargin)
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

            if (~(nargin==2 || nargin==3) || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end

            if (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.targetDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with targetDim() rows.')
            end

            if (nargin<3) accuracy = 1e-6; else accuracy = varargin{2}; end

            uu = EigenMatrix(size(varargin{2},1),size(varargin{2},2),varargin{2});
            result = EigenMatrix();
            calllib('libgismo','gsGeometry_invertPoints',this.objectHandle,uu.ptr(),result.ptr(),accuracy);
            [varargout{1:nargout}] = result.asMatrix();
        end

    end
end
