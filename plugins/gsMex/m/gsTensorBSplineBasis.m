% @file gsTensorBSplineBasis.m
%
%    @brief Matlab wrapper for gsTensorBSplineBasis class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon

classdef gsTensorBSplineBasis < handle

    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
        parDim; % Parametric dimension of the basis
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsTensorBSplineBasis(varargin)
            %gsTensorBSplineBasis - construct a gsTensorBSplineBasis object
            %
            %Usage:
            %  bspb = gsTensorBSplineBasis( file, parDim )
            %  OR
            %  bspb = gsTensorBSplineBasis( knots, parDim )
            %  OR
            %  bspb = gsTensorBSplineBasis( basis, parDim )
            %
            %Input:
            %  parDim: int, parametric dimension of the basis.
            %  file: char, [1 x numChar].
            %    Name of input file from which to read/construct the
            %    gsTensorBSplineBasis.
            %  OR
            %  knots: cell
            %    The cell element i contains the knot vector in the
            %    direction i.
            %  OR
            %  basis: gsTensorBSplineBasis
            %    Copy constructor.
            %
            %Output:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.

            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            if (isa(varargin{1},'uint64'))
                this.objectHandle = varargin{1};
                this.parDim = varargin{2};
            else
                if (~( isa(varargin{1},'char') || isa(varargin{1},'cell') ||...
                        isa(varargin{1},'gsTensorBSplineBasis') ))
                    error(['Input argument no. 1 should be of type ''char'' ',...
                           'or ''gsTHBSplineBasis'' or be a cell of knot vectors.'])
                elseif (isa(varargin{1},'char') && (~exist(varargin{1},'file')))
                    error('File does not exist: %s.',varargin{1})
                end
                if isa(varargin{1},'gsTensorBSplineBasis')
                    var1 = varargin{1}.objectHandle;
                else
                    var1 = varargin{1};
                end
                if (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) ||...
                        ~(floor(varargin{2})==varargin{2}) || ~(varargin{2}>0) )
                    error('Input argument no.2 shoud be a positive integer.')
                else
                    this.parDim = varargin{2};
                end
                this.objectHandle = mex_gsTensorBSplineBasis('constructor', ...
                    class(varargin{1}), var1, this.parDim);
            end
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsTensorBSplineBasis object
            %
            %Usage:
            %  bspb.delete()
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %
            %Output:
            %  (none)
            
            mex_gsTensorBSplineBasis('destructor', this.objectHandle, this.parDim);
        end

        % domainDim - call class method
        function varargout = domainDim(this, varargin)
            %domainDim - dimension of the parameter space of a gsTensorBSplineBasis object
            %
            %Usage:
            %  val = bspb.domainDim()
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %
            %Output:
            %  val: double, [1 x 1].
            %    Dimension of the parameter space of the gsTensorBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('accessor', ...
                this.objectHandle, 'domainDim',  varargin{:}, this.parDim);
        end
        
        % numElements - call class method
        function varargout = numElements(this, varargin)
            %numElements - number of elements of a gsTensorBSplineBasis object
            %
            %Usage:
            %  num = bspb.numElements()
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Number of elements of the gsTensorBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('accessor',...
                this.objectHandle, 'numElements',  varargin{:}, this.parDim);
        end
        
        % size - call class method
        function varargout = size(this, varargin)
            %size - size/number of degrees of freedom of a gsTensorBSplineBasis 
            %   object
            %
            %Usage:
            %  num = bspb.size()
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %
            %Output:
            %  num: double, [1 x 1].
            %    Size of the gsTensorBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('accessor', ...
                this.objectHandle, 'size',  varargin{:}, this.parDim);
        end

        % support - call class method
        function varargout = support(this, varargin)
            %support - support of a gsTensorBSplineBasis object
            %
            %Usage:
            %  supp = bspb.support()
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %
            %Output:
            %  supp: double,  [d x 2].
            %    Support of the gsTensorBSplineBasis, ordered like 
            %      [u1_min, u1_max; u2_min, u2_max; ..., ud_min, ud_max]
            %    where d is the parametric dimension of the 
            %    gsTensorBSplineBasis.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('accessor', ...
                this.objectHandle, 'support',  varargin{:}, this.parDim);
        end

        % degree - call class method
        function [varargout] = degree(this, varargin)
            %degree - the degree for a specified direction of a 
            %   gsTensorBSplineBasis object
            %
            %Usage:
            %  deg = bspb.degree( dir )
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  dir: int, [1 x 1].
            %    Direction of space for which we want to know the degree of 
            %    the gsTensorBSplineBasis.
            %
            %Output:
            %  deg: double, [1 x 1].
            %    Degree of the gsTensorBSplineBasis object in direction dir. 
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}<1 || ...
                    varargin{1}>this.domainDim())
                error('Input argument must be a non negative integer less than %d.',...
                    this.domainDim())
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('degree', ...
                this.objectHandle, varargin{:}, this.parDim);
        end
        
        % eval - call class method
        function [varargout] = eval(this, varargin)
            %eval - evaluate a gsTensorBSplineBasis object
            %
            %Usage:
            %  val = bspb.eval( pts )
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  pts: double, [d x numPts] or [numPts x d].
            %    Points in which to evaluate the gsTensorBSplineBasis.
            %
            %Output:
            %  val: double, [numFun x numPts].
            %    Value of all active functions in each of the specified
            %    points.
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}))
                error('Input argument no. 1 must be numeric, 2-dimensional.')
            elseif (~isequal(size(varargin{1},2),this.domainDim()) && ...
                    ~isequal(size(varargin{1},1),this.domainDim()))
                error(['Input argument no. 1 must be numeric, 2-dimensional ',...
                    'with one dimension equal to parametric dim.'])
            elseif (isequal(size(varargin{1},2),this.domainDim()))
                varargin{1} = varargin{1}.';
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('eval', ...
                this.objectHandle, varargin{:}, this.parDim);
        end

        % evalSingle - call class method
        function [varargout] = evalSingle(this, varargin)
            %evalSingle - evaluate a single function in a gsTensorBSplineBasis object
            %
            %Usage:
            %  valSingle = bspb.evalSingle( fun, pts )
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
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
            elseif (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) ||...
                    ~isequal(size(varargin{2},1),this.domainDim()))
                error('Input argument no. 2 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('evalSingle', ...
                this.objectHandle, varargin{:}, this.parDim);
        end
        
        % save - call class method
        function [varargout] = save(this, varargin)
            %save - save a gsTensorBSplineBasis object as xml object
            %
            %Usage:
            %  bspb.save(filename);
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
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
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('save', ...
                this.objectHandle, varargin{:}, this.parDim);
        end
        
        % knots - call class method
        function [varargout] = knots(this, varargin)
            %knots - returns the knot vector of a gsTensorBSplineBasis object
            %   on the specified direction
            %
            %Usage:
            %  knt = bspb.knots( dir )
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  dir: int, [1 x 1].
            %    index of the direction of space to consider. 
            %
            %Output:
            %  knt: double, [1 x numKnots].
            %    Knot vector in the direction dir. 
            
            if (nargin~=2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                    ~(mod(varargin{1},1)==0) || varargin{1}<1 || ...
                    varargin{1}>this.domainDim)
                error('Input argument no. 1 must be a strictly positive integer smaller than %d.',...
                    this.domainDim)
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('knots', ...
                this.objectHandle, varargin{:}, this.parDim);
        end
        
        % active - call class method
        function [varargout] = active(this, varargin)
            %active - active functions of a gsTensorBSplineBasis object
            %
            %Usage:
            %  act = bspb.active( pts )
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
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
                    ~isequal(size(varargin{1},1),this.domainDim()))
                error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis('active', ...
                this.objectHandle, varargin{:}, this.parDim);
        end

        % uniformRefine - call class method
        function [varargout] = uniformRefine(this, varargin)
            %uniformRefine - proceed to a uniformRefinement of a gsTensorBSplineBasis object
            %
            %Usage:
            %  bspb.uniformRefine(numKnots, mult)
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  numKnots: int, number of new knots inserted on each knot span.
            %  mult: int, each knot is inserted with a multiplicity mult.
            %
            %Output:
            %  None - changes done inline.

            if (nargin~=3 || nargout>0)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
                ~(floor(varargin{1})==varargin{1}) || ~(varargin{1}>0) ||...
                ~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || ...
                ~(floor(varargin{2})==varargin{2}) || ~(varargin{2}>0))
                error('Input argument no. 1 and 2 must stricly positive integers.')
            end
            mex_gsTensorBSplineBasis('uniformRefine', this.objectHandle, ...
                varargin{:}, this.parDim);
        end

        % uniformRefine_withCoefs - call class method
        function [varargout] = uniformRefine_withCoefs(this, varargin)
            %uniformRefine_withCoefs - Refine the basis given by a gsTensorBSplineBasis 
            % object uniformly. The function simultainously updates the vector coefs,
            % representing a function in the bases, such that its new version 
            % represents the same function.
            %
            %Usage:
            %  coefs_ref = bspb.uniformRefine_withCoefs(coefs, numKnots, mult)
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  coefs: the coefficients representing a function in the bases 
            %    before the basis is refined.
            %  numKnots: int, number of new knots inserted on each knot span.
            %  mult: int, each knot is inserted with a multiplicity mult.
            %
            %Output:
            %  coefs_ref: the coefficients representing a function in the
            %     bases after the basis has been refined.

            if (nargin~=4 || nargout~=1)
                error('Invalid number of input and/or output arguments.')
            end
            if (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) ||...
                ~(floor(varargin{2})==varargin{2}) || ~(varargin{2}>0) ||...
                ~isa(varargin{3},'numeric') || ~isscalar(varargin{3}) ||...
                ~(floor(varargin{3})==varargin{3}) || ~(varargin{3}>0) ||...
                ~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}))
                error(['Input argument no. 1 should be a 2d array of double, ',...
                    'inputs no. 2 and 3 must be stricly positive integers.'])
            end
            [varargout{1:nargout}] = mex_gsTensorBSplineBasis(...
                'uniformRefine_withCoefs', this.objectHandle,...
                 varargin{:}, this.parDim);
        end

        % refineElements - call class method
        function [varargout] = refineElements(this, varargin)
            %refineElements - proceed to a refinement of given elements of 
            % a gsTensorBSplineBasis object
            %
            %Usage:
            %  bspb.refineElements(boxes)
            %
            %Input:
            %  bspb: gsTensorBSplineBasis, [1 x 1].
            %    The gsTensorBSplineBasis object.
            %  boxes: int array of dimension N(2d+1) where N is the number
            %    of boxes to refine (sets of cells), and d is the parametric
            %    dimension. 
            %
            %Output:
            %  None - changes done inline.

            if (nargin~=2 || nargout>0)
                error('Invalid number of input and/or output arguments.')
            end
            if ( ~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ...
                    ~(min(size(varargin{1}))==1) )
                error('Input argument must be an array of integers.')
            end
            mex_gsTensorBSplineBasis('refineElements', this.objectHandle, ...
                varargin{:}, this.parDim);
        end

%         VIRTUAL METHOD IN G+SMO
%         % refineElements_withCoefs - call class method
%         function [varargout] = refineElements_withCoefs(this, varargin)
%             %refineElements_withCoefs - Refine the given elements of the 
%             % basis given by a gsTensorBSplineBasis object. The function 
%             % simultainously updates the vector coefs, representing a 
%             % function in the bases, such that its new version represents 
%             % the same function.
%             %
%             %Usage:
%             %  coefs_ref = bspb.refineElements_withCoefs(coefs, boxes)
%             %
%             %Input:
%             %  bspb: gsTensorBSplineBasis, [1 x 1].
%             %    The gsTensorBSplineBasis object.
%             %  coefs: the coefficients representing a function in the bases 
%             %    before the basis is refined.
%             %  boxes: int array of dimension N(2d+1) where N is the number
%             %    of boxes to refine (sets of cells), and d is the parametric
%             %    dimension. 
%             %
%             %Output:
%             %  coefs_ref: the coefficients representing a function in the
%             %    bases after the basis has been refined.
% 
%             if (nargin~=3 || nargout~=1)
%                 error('Invalid number of input and/or output arguments.')
%             end
%             if ( ~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) ||...
%                     ~(min(size(varargin{2}))==1) ||...
%                 ~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) )
%                 error(['Input argument no. 1 should be a 2d array of double, ',...
%                     'inputs no. 2 must be an array of integers.'])
%             end
%             [varargout{1:nargout}] = mex_gsTensorBSplineBasis(...
%                 'refineElements_withCoefs', this.objectHandle,...
%                 varargin{:}, this.parDim);
%         end
        
    end
end
