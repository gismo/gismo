% @file gsTensorBSpline.m
%
%    @brief Matlab wrapper for gsTensorBSpline class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon

classdef gsTensorBSpline < handle

    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
        paramDim; % Parametric dimension
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsTensorBSpline(varargin)
            %gsTensorBSpline - construct a gsTensorBSpline object
            %
            %Usage:
            %  bsp = gsTensorBSpline( file, paramDim )
            %  OR
            %  bsp = gsTensorBSpline( basis, coefs, paramDim )
            %  OR
            %  bsp = gsTensorBSpline( geom, paramDim )
            %
            %Input:
            %  paramDim: int, parametric dimension of the basis.
            %  file: char, [1 x numChar].
            %    Name of input file from which to read/construct the
            %    gsTensorBSpline.
            %  OR
            %  bbasis: gsTensorBSplineBasis
            %    Bspline basis from which the geometry is built
            %  coefs: array of double of size [numCoefs x geoDim],
            %    where numCoefs is the number of coefficients (total
            %    number of active basis functions) and geoDim is the
            %    dimension of the physical space.
            %  OR
            %  geom: gsTensorBSpline
            %    Copy constructor.
            %
            %Output:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            
            if (nargin>3 || nargin<2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            elseif (nargin==2)
                if isa(varargin{1},'uint64')
                    this.objectHandle = varargin{1};
                    this.paramDim = varargin{2};
                else
                    if (~(isa(varargin{1},'char') || isa(varargin{1},'gsTensorBSpline') ))
                        error(['First input arguments should be of type ''char'',', ...
                               'or a gsTensorBSpline, or a gsTensorBSplineBasis and a ',...
                               '2d-array of double.'])
                    elseif (isa(varargin{1},'char') && ~exist(varargin{1},'file'))
                        error('File does not exist: %s.',varargin{1})
                    end
                    if isa(varargin{1}, 'gsTensorBSpline')
                        var1 = varargin{1}.objectHandle;
                    else
                        var1 = varargin{1};
                    end
                    if (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) ||...
                        ~(floor(varargin{2})==varargin{2}) || ~(varargin{2}>0) )
                        error('Input argument no.2 shoud be a positive integer.')
                    else
                        this.paramDim = varargin{2};
                    end
                    this.objectHandle = mex_gsTensorBSpline('constructor', ...
                        class(varargin{1}), var1, this.paramDim);
                end
            elseif (nargin==3)
                var2 = varargin{2};
                if (~(isa(varargin{1},'gsTensorBSplineBasis') && isa(var2,'double') && ismatrix(var2)))
                    error(['Input arguments should be of type ''char'', ',...
                        'or a gsTensorBSplineBasis and a 2d-array of double.'])
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
                this.objectHandle = mex_gsTensorBSpline('constructor', class(varargin{1}),...
                    class(varargin{2}), var1, var2, this.paramDim);
            end
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsTensorBSpline object
            %
            %Usage:
            %  bsp.delete()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  (none)
            
            mex_gsTensorBSpline('destructor', this.objectHandle, this.paramDim);
        end

        % parDim - call class method
        function varargout = parDim(this, varargin)
            %parDim - dimension of the parameter space of a gsTensorBSpline object
            %
            %Usage:
            %  val = bsp.parDim()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the parameter space of the gsTensorBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('accessor', ...
                this.objectHandle, 'parDim',  varargin{:}, this.paramDim);
        end
        
        % geoDim - call class method
        function varargout = geoDim(this, varargin)
            %geoDim - dimension of the physical space of a gsTensorBSpline object
            %
            %Usage:
            %  val = bsp.geoDim()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the physical space of the gsTensorBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('accessor', ...
                this.objectHandle, 'geoDim',  varargin{:}, this.paramDim);
        end

        % size - call class method
        function varargout = size(this, varargin)
            %size - size of a gsTensorBSpline object
            %
            %Usage:
            %  num = bsp.size()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the gsTensorBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('accessor', ...
                this.objectHandle, 'size',  varargin{:}, this.paramDim);
        end
        
        % support - call class method
        function varargout = support(this, varargin)
            %support - support of a gsTensorBSpline object
            %
            %Usage:
            %  supp = bsp.support()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  supp: double, [d x 2].
            %    Support of the gsTensorBSpline, ordered like 
            %      [u1_min, u1_max; u2_min, u2_max; ..., ud_min, ud_max]
            %    where d is the parametric dimension of the 
            %    gsTensorBSpline.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('accessor', ...
                this.objectHandle, 'support',  varargin{:}, this.paramDim);
        end
        
        % basis - call class method
        function [varargout] = basis(this, varargin)
            %basis - returns the gsTensorBSplineBasis object linked to a 
            % gsTensorBSpline object
            %
            %Usage:
            %  act = bsp.basis()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  basis: gsTensorBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            basis_ptr = mex_gsTensorBSpline('accessor', this.objectHandle, ...
                'basis', varargin{:}, this.paramDim);
            [varargout{1:nargout}] = gsTensorBSplineBasis(basis_ptr, this.paramDim);
        end
        
        % coefs - call class method
        function varargout = coefs(this, varargin)
            %coefs - coefficients/control points of a gsTensorBSpline object
            %
            %Usage:
            %  supp = bsp.coefs()
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %
            %Output:
            %  cc: array of double. Control points of the gsTensorBSpline.
            %    of size [numCoefs x geoDim], where numCoefs is the 
            %    number of coefficients and geoDim is the dimension 
            %    of the physical space.

            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('accessor', ...
                this.objectHandle, 'coefs', varargin{:}, this.paramDim);
        end

        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsTensorBSpline object
            %
            %Usage:
            %  val = bsp.eval( pts )
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTensorBSpline.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) ||...
                    ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('eval', ...
                this.objectHandle, varargin{:}, this.paramDim);
        end
        
        % jacobian - call class method
        function [varargout] = jacobian(this, varargin)
            %jacobian - evaluate the jacobian of a gsTensorBSpline object
            %
            %Usage:
            %  val = bsp.jacobian( pts )
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTensorBSpline.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of the jacobian of all active functions in each of
            %    the specified points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) ||...
                    ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('jacobian', ...
                this.objectHandle, varargin{:}, this.paramDim);
        end

        % hess - call class method
        function [varargout] = hess(this, varargin)
            %deriv - evaluate the hessian of a gsTensorBSpline object in one
            %direction
            %
            %Usage:
            %  val = bsp.hess( pts, dir )
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the gsTensorBSpline.
            %  dir: int
            %    Direction of space on which to compute the hessian.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of the hessian matrix in direction dir of all active 
            %    functions in each of the specified points.
            
            if (nargin~=3 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            elseif (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) ||...
                    ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with %d rows.',...
                    this.parDim())
            elseif (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || ...
                    ~(mod(varargin{2},1)==0) || varargin{2}<1 || varargin{2}>this.parDim())
                error('Input argument no. 2 must be a non negative integer smaller than %d.',...
                    this.parDim())
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('hess', ...
                this.objectHandle, varargin{:}, this.paramDim);
        end
        
        % active - call class method
        function [varargout] = active(this, varargin)
            %active - active functions of a gsTensorBSpline object
            %
            %Usage:
            %  act = bsp.active( pts )
            %
            %Input:
            %  bsp: gsTensorBSpline, [1 x 1].
            %    The gsTensorBSpline object.
            %  pts: double, [d x numPts].
            %    Points in which to evaluate the function.
            %
            %Output:
            %  act: double, [numFun x numPts].
            %    Index of active functions in each of the specified points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ...
                    ~isequal(size(varargin{1},1),this.parDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSpline('active', ...
                this.objectHandle, varargin{:}, this.paramDim);
        end
        
        % save - call class method
        function [varargout] = save(this, varargin)
            %save - save a gsTensorBSpline object as xml object
            %
            %Usage:
            %  bsp.save(filename);
            %
            %Input:
            %  bsp: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
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
            [varargout{1:nargout}] = mex_gsTensorBSpline('save', ...
                this.objectHandle, varargin{:}, this.paramDim);
        end

    end
end
