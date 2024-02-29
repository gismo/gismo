% @file gsKnotVector.m
%
%    @brief Matlab wrapper for gsKnotVector class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): O. Chanon, P. Noertoft

classdef gsKnotVector < handle
% try: Custom Display Interface ? override disp ?
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsKnotVector(varargin)
            %gsKnotVector - construct a gsKnotVector object
            %
            %Usage:
            %  gsKv = gsKnotVector( deg,kv )
            %
            %Input:
            %  deg: int, [1x1]
            %    Degree of the knot vector
            %  kv: real, [numKnots x 1].
            %    Row vector representing the knots
            %
            %Output:
            %  gsKv: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            % HV: WHAT DOES THIS ONE DO??
            if (isa(varargin{1},'uint64'))
                this.objectHandle = varargin{1};
            else
                if (~(isa(varargin{1},'double')))
                    error('Input argument no. 2 should be of type ''double''.')
                end
                this.objectHandle = mex_gsKnotVector('constructor', class(varargin{1}), varargin{:});
            end
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            %delete - delete a gsKnotVector object
            %
            %Usage:
            %  thb.delete()
            %
            %Input:
            %  thb: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            %
            %Output:
            %  (none)
        
            mex_gsKnotVector('destructor', this.objectHandle);
        end

         % This function returns the address of the C++ pointer
         function varargout = ptr(this, varargin)
             if (nargin~=1 || nargout>1)
                 error('Invalid number of input and/or output arguments.')
             end
             varargout{1} = this.objectHandle;
         end

        function disp(this, var_name)
            a = mex_gsKnotVector('get',this.objectHandle,[]);
            %with pointer address:
            %fmt=['gsknotVector(0x%s):' repmat(' %1.0f',1,numel(a)) '\n'];
            %fprintf(fmt,dec2hex(this.objectHandle),a);
            fmt=['gsknotVector {' repmat(' %1.0f',1,numel(a)) ' }\n'];
            fprintf(fmt,a);
        end

        % % dim - call class method
        % function varargout = dim(this, varargin)
        %     %dim - dimension of the parameter space of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  val = thb.dim()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  val: double, [1 x 1].
        %     %    Dimension of the parameter space of the gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'dim',  varargin{:});
        % end
        
        % % numElements - call class method
        % function varargout = numElements(this, varargin)
        %     %numElements - number of elements of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  num = thb.numElements()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  num: double, [1 x 1].
        %     %    Number of elements of the gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'numElements',  varargin{:});
        % end
        
        % % support - call class method
        % function varargout = support(this, varargin)
        %     %support - support of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  supp = thb.support()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  supp: double, [1 x 2*d].
        %     %    Support of the gsKnotVector, ordered like
        %     %      [u1_min, ..., ud_min, u1_max, ..., ud_max]
        %     %    where d is the parametric dimennsion of the
        %     %    gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'support',  varargin{:});
        % end

        % % size - call class method
        % function varargout = size(this, varargin)
        %     %size - size of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  num = thb.size()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  num: double, [1 x 1].
        %     %    Size of the gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'size',  varargin{:});
        % end

        % % treeSize - call class method
        % function varargout = treeSize(this, varargin)
        %     %treeSize - size of the tree of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  num = thb.treeSize()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  num: double, [1 x 1].
        %     %    Size of the tree of the gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'treeSize',  varargin{:});
        % end
        
        % % treeLeafSize - call class method
        % function varargout = treeLeafSize(this, varargin)
        %     %treeLeafSize - size of the leaf in the tree of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  num = thb.treeLeafSize()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  num: double, [1 x 1].
        %     %    Size of the leaf in the tree of the gsKnotVector.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'treeLeafSize',  varargin{:});
        % end

        % % maxLevel - call class method
        % function varargout = maxLevel(this, varargin)
        %     %maxLevel - maximum level of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  lev = thb.maxLevel()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  lev: double, [1 x 1].
        %     %    Maximum level present in the hierarchy of the
        %     %    gsKnotVector object.
            
        %     if (nargin~=1 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('accessor', this.objectHandle, 'maxLevel',  varargin{:});
        % end
        
        % % treePrintLeaves - call class method
        % function varargout = treePrintLeaves(this, varargin)
        %     %treePrintLeaves - print the leaves in the tree of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  thb.treePrintLeaves()
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %  (none - outputs to the screen).
            
        %     if (nargin~=1 || nargout>0)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('treePrintLeaves', this.objectHandle, varargin{:});
        % end

        % degree - call class method
        function [varargout] = degree(this, varargin)
            %degree - the degree for a specified direction of a 
            %   gsKnotVector object
            %
            %Usage:
            %  deg = thb.degree( dir )
            %
            %Input:
            %  thb: gsKnotVector, [1 x 1].
            %    The gsKnotVector object.
            %  dir: int, [1 x 1].
            %    Direction of space for which we want to know the degree of 
            %    the gsKnotVector.
            %
            %Output:
            %  deg: double, [1 x 1].
            %    Degree of the gsKnotVector object in direction dir.
            
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = mex_gsKnotVector('degree', this.objectHandle, varargin{:});
        end
        
        % % eval - call class method
        % function [varargout] = eval(this, varargin)
        %     %eval - evaluate a gsKnotVector object
        %     %
        %     %Usage:
        %     %  val = thb.eval( pts )
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %  pts: double, [d x numPts].
        %     %    Points in which to evaluate the gsKnotVector.
        %     %
        %     %Output:
        %     %  val: double, [numFun x numPts].
        %     %    Value of all active functions in each of the specified
        %     %    points.
            
        %     if (nargin~=2 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
        %         error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('eval', this.objectHandle, varargin{:});
        % end

        % % evalSingle - call class method
        % function [varargout] = evalSingle(this, varargin)
        %     %evalSingle - evaluate a single function in a gsKnotVector object
        %     %
        %     %Usage:
        %     %  valSingle = thb.evalSingle( fun, pts )
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %  fun: double, [1 x 1].
        %     %    Index of function to evaluate.
        %     %  pts: double, [d x numPts].
        %     %    Points in which to evaluate the function.
        %     %
        %     %Output:
        %     %  val: double, [1 x numPts].
        %     %    Value of the specified function in each of the specified
        %     %    points.
            
        %     if (nargin~=3 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ...
        %             ~(mod(varargin{1},1)==0) || varargin{1}<1)
        %         error('Input argument no. 1 must be an strictly positive integer.')
        %     elseif (~isa(varargin{2},'numeric') || ~ismatrix(varargin{2}) || ~isequal(size(varargin{2},1),this.dim()))
        %         error('Input argument no. 2 must be numeric, 2-dimensional, and with d rows.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('evalSingle', this.objectHandle, varargin{:});
        % end
        
        % % save - call class method
        % function [varargout] = save(this, varargin)
        %     %save - save a gsKnotVector object as xml object
        %     %
        %     %Usage:
        %     %  thb.save();
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %
        %     %Output:
        %     %   (none - saved into an xml file)
            
        %     if (nargin~=2)
        %         error('Invalid number of input arguments.')
        %     end
        %     if (~(isa(varargin{1},'char')))
        %         error('Input argument no. 1 should be of type ''char''.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('save', this.objectHandle, varargin{:});
        % end
        
        % % knots - call class method
        % function [varargout] = knots(this, varargin)
        %     %knots - returns the knot vector of a gsKnotVector object
        %     %   of the specified level on the specified direction
        %     %
        %     %Usage:
        %     %  knt = thb.knots( lev, dir )
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %  lev: int, [1 x 1].
        %     %    Index of the level of the hierarchy to consider.
        %     %  dir: int, [1 x 1].
        %     %    index of the direction of space to consider.
        %     %
        %     %Output:
        %     %  knt: double, [1 x numKnots].
        %     %    Knot vector corresponding to level lev in the direction
        %     %    dir.
            
        %     if (nargin~=3 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     if (~isa(varargin{1},'numeric') || ~isscalar(varargin{1}) || ~(mod(varargin{1},1)==0) || varargin{1}<1)
        %         error('Input argument no. 1 must be a strictly positive integer.')
        %     elseif (~isa(varargin{2},'numeric') || ~isscalar(varargin{2}) || ...
        %             ~(mod(varargin{2},1)==0) || varargin{2}<1 || varargin{2}>this.dim())
        %         error('Input argument no. 2 must be a strictly positive integer smaller than %d.', this.dim())
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('knots', this.objectHandle, varargin{:});
        % end
        
        % % active - call class method
        % function [varargout] = active(this, varargin)
        %     %active - active functions of a gsKnotVector object
        %     %
        %     %Usage:
        %     %  act = thb.active( pts )
        %     %
        %     %Input:
        %     %  thb: gsKnotVector, [1 x 1].
        %     %    The gsKnotVector object.
        %     %  pts: double, [d x numPts].
        %     %    Points in which to evaluate the function.
        %     %
        %     %Output:
        %     %  act: double, [numFun x numPts].
        %     %    Index of active functions in each of the specified points.
            
        %     if (nargin~=2 || nargout>1)
        %         error('Invalid number of input and/or output arguments.')
        %     end
        %     if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}) || ~isequal(size(varargin{1},1),this.dim()))
        %         error('Input argument no. 1 must be numeric, 2-dimensional, and with d rows.')
        %     end
        %     [varargout{1:nargout}] = mex_gsKnotVector('active', this.objectHandle, varargin{:});
        % end

    end
end
