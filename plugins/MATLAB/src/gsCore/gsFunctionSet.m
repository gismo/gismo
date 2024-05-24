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
%    Author(s): H.M. Verhelst

classdef gsFunctionSet < handle

    % properties (SetAccess = private, Hidden = true)
    %     objectHandle; % Handle to the underlying C++ class instance
    % end

    methods(Access = public)

         % This function returns the address of the C++ pointer
         function varargout = ptr(this, varargin)
             if (nargin~=1 || nargout>1)
                 error('Invalid number of input and/or output arguments.')
             end
             varargout{1} = this.objectHandle;
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

            calllib('libgismo','gsFunctionSet_delete',this.objectHandle);
        end

        function disp(this, var_name)
            calllib('libgismo','gsFunctionSet_print',this.objectHandle);
        end

        % domainDim - call class method
        function varargout = domainDim(this, varargin)
            %dim - dimension of the parameter space of a gsFunctionSet object
            %
            %Usage:
            %  val = f.domainDim()
            %
            %Input:
            %  f: gsFunctionSet, [1 x 1].
            %    The gsFunctionSet object.
            %
            %Output:
            %  val: int, [1 x 1].
            %    Dimension of the parameter space of the gsBasis.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsFunctionSet_domainDim',this.objectHandle);
        end

        % targetDim - call class method
        function varargout = targetDim(this, varargin)
            %dim - dimension of the target space of a gsFunctionSet object
            %
            %Usage:
            %  val = f.targetDim()
            %
            %Input:
            %  f: gsFunctionSet, [1 x 1].
            %    The gsFunctionSet object.
            %
            %Output:
            %  val: int, [1 x 1].
            %    Dimension of the target space of the gsBasis.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsFunctionSet_targetDim',this.objectHandle);
        end

        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsFunctionSet object
            %
            %Usage:
            %  val = thb.eval( pts )
            %
            %Input:
            %  thb: gsFunctionSet, [1 x 1].
            %    The gsFunctionSet object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsFunctionSet.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.

            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.domainDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{1},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrix();
            calllib('libgismo','gsFunctionSet_eval_into',this.objectHandle,uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

        % deriv - call class method
        function [varargout] = deriv(this, varargin)
            %deriv - take the derivative of a gsFunctionSet object
            %
            %Usage:
            %  val = thb.deriv( pts )
            %
            %Input:
            %  thb: gsFunctionSet, [1 x 1].
            %    The gsFunctionSet object.
            %  pts: double, [d x numPts].
            %    Points in which to take the derivative of the gsFunctionSet.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.

            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.domainDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{1},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrix();
            calllib('libgismo','gsFunctionSet_deriv_into',this.objectHandle,uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end

        % deriv2 - call class method
        function [varargout] = deriv2(this, varargin)
            %deriv2 - take the second derivative of a gsFunctionSet object
            %
            %Usage:
            %  val = thb.deriv2( pts )
            %
            %Input:
            %  thb: gsFunctionSet, [1 x 1].
            %    The gsFunctionSet object.
            %  pts: double, [d x numPts].
            %    Points in which to take the second derivative of the gsFunctionSet.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.

            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.domainDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end

            assert(size(varargin{1},1)==this.domainDim(),"Domain dimension should be equal to the number of rows of the points");
            uu = EigenMatrix(size(varargin{1},1),size(varargin{1},2),varargin{1});
            result = EigenMatrix();
            calllib('libgismo','gsFunctionSet_deriv2_into',this.objectHandle,uu.ptr(),result.ptr());
            [varargout{1:nargout}] = result.asMatrix();
        end
    end
end