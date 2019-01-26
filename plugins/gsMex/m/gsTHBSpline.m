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
        paramDim; % Parametric dimension
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsTHBSpline(varargin)
            %gsTHBSpline - construct a gsTHBSpline object
            %
            %Usage:
            %  thb = gsTHBSpline( file, paramDim )
            %  OR
            %  thb = gsTHBSpline( thbbasis, coefs, paramDim )
            %  OR
            %  thb = gsTHBSpline( thbgeom, paramDim )
            %
            %Input:
            %  paramDim: int, parametric dimension.
            %  file: char, [1 x numChar].
            %    Name of input file from which to read/construct the
            %    gsTHBSpline.
            %  OR
            %  thbbasis: gsTHBSplineBasis
            %    THB spline basis from which the geometry is built
            %  coefs: array of double of size [numCoefs x geoDim],
            %    where numCoefs is the number of coefficients (total
            %    number of active basis functions) and geoDim is the
            %    dimension of the physical space.
            %  OR
            %  thbgeom: gsTHBSpline
            %    Copy constructor.
            %
            %Output:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            
            if (nargin>3 || nargin<2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            elseif (nargin==2)
                if isa(varargin{1},'uint64')
                    this.objectHandle = varargin{1};
                    this.paramDim = varargin{2};
                else
                    if (~(isa(varargin{1},'char') || isa(varargin{1},'gsTHBSpline') ))
                        error(['First input arguments should be of type ''char'',', ...
                               'or a gsTHBSpline, or a gsTHBSplineBasis and a ',...
                               '2d-array of double.'])
                    elseif (isa(varargin{1},'char') && ~exist(varargin{1},'file'))
                        error('File does not exist: %s.',varargin{1})
                    end
                    if isa(varargin{1}, 'gsTHBSpline')
                        var1 = varargin{1}.objectHandle;
                    else
                        var1 = varargin{1};
                    end
                    if (~isa(varargin{2},'numeric') || ~isscalar(varargin{2})...
                            || ~(floor(varargin{2})==varargin{2}) || ~(varargin{2}>0) )
                        error('Last input argument shoud be a positive integer.')
                    else
                        this.paramDim = varargin{2};
                    end
                    this.objectHandle = mex_gsTHBSpline('constructor', ...
                        class(varargin{1}), var1, this.paramDim);
                end
            elseif (nargin==3)
                var2 = varargin{2};
                if (~(isa(varargin{1},'gsTHBSplineBasis') && isa(var2,'double') && ismatrix(var2)))
                    error(['First input arguments should be of type ''char'', '...
                           'or a gsTHBSplineBasis and a 2d-array of double.'])
                % elseif (size(var2,1)~=varargin{1}. TODO!!! number of dof!)
                %    error('Wrong coefficient dimension with respect to the basis.')
                end
                if (~isa(varargin{3},'numeric') || ~isscalar(varargin{3}) || ...
                        ~(floor(varargin{3})==varargin{3}) || ~(varargin{3}>0) )
                    error('Last input argument shoud be a positive integer.')
                else
                    this.paramDim = varargin{3};
                end
                var1 = varargin{1}.objectHandle;
                this.objectHandle = mex_gsTHBSpline('constructor', class(varargin{1}),...
                    class(varargin{2}), var1, var2, this.paramDim);
            end
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
            
            mex_gsTHBSpline('destructor', this.objectHandle, this.paramDim);
        end

        % parDim - call class method
        function varargout = parDim(this, varargin)
            %parDim - dimension of the parameter space of a gsTHBSpline object
            %
            %Usage:
            %  val = thb.parDim()
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
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'parDim',  varargin{:}, this.paramDim);
        end
        
        % geoDim - call class method
        function varargout = geoDim(this, varargin)
            %geoDim - dimension of the physical space of a gsTHBSpline object
            %
            %Usage:
            %  val = thb.geoDim()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the physical space of the gsTHBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'geoDim',  varargin{:}, this.paramDim);
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
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'size',  varargin{:}, this.paramDim);
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
            %  supp: double, [d x 2].
            %    Support of the gsTHBSpline, ordered like 
            %      [u1_min, u1_max; u2_min, u2_max; ..., ud_min, ud_max]
            %    where d is the parametric dimension of the 
            %    gsTHBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'support',  varargin{:}, this.paramDim);
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
            basis_ptr = mex_gsTHBSpline('accessor', this.objectHandle, 'basis', varargin{:}, this.paramDim);
            [varargout{1:nargout}] = gsTHBSplineBasis(basis_ptr, this.paramDim);
        end
        
        % coefs - call class method
        function varargout = coefs(this, varargin)
            %coefs - coefficients/control points of a gsTHBSpline object
            %
            %Usage:
            %  supp = thb.coefs()
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %
            %Output:
            %  cc: array of double. Control points of the gsTHBSpline.
            %    of size [numCoefs x geoDim], where numCoefs is the 
            %    number of coefficients (total number of active basis 
            %    functions) and geoDim is the dimension of the physical space.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('accessor', this.objectHandle, 'coefs', varargin{:}, this.paramDim);
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
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.paramDim))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('eval', this.objectHandle, varargin{:}, this.paramDim);
        end
        
        % jacobian - call class method
        function [varargout] = jacobian(this, varargin)
            %jacobian - evaluate the jacobian of a gsTHBSpline object
            %
            %Usage:
            %  val = thb.jacobian( pts )
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTHBSpline.
            %
            %Output:
            %  val: double, [geoDim x parDim x numPts].
            %    Value of the jacobian of all active functions in each of
            %    the specified points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('jacobian', this.objectHandle, varargin{:}, this.paramDim);
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
            %  val: double, [geoDim x parDim x parDim x numPts].
            %    Value of the hessian matrix in direction dir of all active 
            %    functions in each of the specified points.
            
            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            elseif (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with %d rows.', this.parDim())
            elseif (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || ...
                    ~(mod(varargin{2},1)==0) || varargin{2}<1 || varargin{2}>this.parDim())
                error('Input argument no. 2 must be a non negative integer smaller than %d.', this.parDim())
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('hess', this.objectHandle, varargin{:}, this.paramDim);
        end

        % sliceCoefs - call class method
        function [varargout] = sliceCoefs(this, varargin)
            %sliceCoefs - cefficients corresponding to an isoparametric slice
            % of this gsTHBSpline object.
            %
            %Usage:
            %  slCoefs = thb.sliceCoefs(dir_fixed, par)
            %
            %Input:
            %  thb: gsTHBSpline, [1 x 1].
            %    The gsTHBSpline object.
            %  dir_fixed: int, direction fixed for slicing.
            %  par: double, parameter fixed for slicing.
            %
            %Output:
            %  slCoefs: array of double [numCoefs x parDim]
            %    The coefficients corresponding to the slice.

            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            elseif (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}<1 || varargin{1}>this.parDim)
                error('Input argument no. 1 must be a non negative integer smaller than %d.', this.parDim)
            elseif (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || numel(varargin{1})~=1)
                error('Input argument no.2 must be a scalar.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('sliceCoefs', this.objectHandle, varargin{:}, this.paramDim);
        end
        
        % save - call class method
        function [varargout] = save(this, varargin)
            %save - save a gsTHBSpline object as xml object
            %
            %Usage:
            %  thb.save(filename);
            %
            %Input:
            %  thb: gsTHBSplineBasis, [1 x 1].
            %    The gsTHBSplineBasis object.
            %  filename: string, name of the xml output file without
            %    extension
            % 
            %Output:
            %   (none - saved into an xml file)
            
            if (nargin~=2)
                error('Invalid number of input arguments.')
            end
            if (~(isa(varargin{1},'char')))
                error('Input argument no. 1 should be of type ''char''.')
            end
            [varargout{1:nargout}] = mex_gsTHBSpline('save', this.objectHandle, varargin{:}, this.paramDim);
        end

    end
end
