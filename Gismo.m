classdef Gismo
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    methods (Static)

        function r = rows(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (class(varargin{1})=='EigenMatrix')
                r = calllib('libgismo','gsMatrix_rows',varargin{1}.ptr);
            else
                error('Data type not understood')
            end
        end

        function c = cols(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (class(varargin{1})=='EigenMatrix')
                c = calllib('libgismo','gsMatrix_cols',varargin{1}.ptr);
            else
                error('Data type not understood')
            end
        end

        function r = data(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (class(varargin{1})=='EigenMatrix')
                r = calllib('libgismo','gsMatrix_data',varargin{1}.ptr);
            else
                error('Data type not understood')
            end
        end

        function ddim = domainDim(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (class(varargin{1})=='Geometry')
                ddim = calllib('libgismo','domainDim',varargin{:});
            else
                error('Data type not understood')
            end
        end

        function tdim = targetDim(varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (class(varargin{1})=='Geometry')
                tdim = calllib('libgismo','targetDim',{varargin{1}.ptr});
            else
                error('Data type not understood')
            end
        end

        function ev = eval(obj,pts)
            if (class(obj)=='Geometry')
                result = EigenMatrix();
                points = EigenMatrix(size(pts,1),size(pts,2),libpointer('doublePtr',pts));
                class(obj)
                class(points)
                class(result)
                calllib('libgismo','eval_into',{obj,points,result});
                ev = result;
            else
                error('Data type not understood')
            end
        end

    end
end