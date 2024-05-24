classdef EigenMatrixInt < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties (SetAccess = public, Hidden = true)
        objectHandle
    end

    methods
        function this = EigenMatrixInt(varargin)
            if (nargin==0)
                this.objectHandle = calllib('libgismo','gsMatrixInt_create');
            elseif(nargin==1)
                if (~isa(varargin{1},'numeric') || ~ismatrix(varargin{1}))
                    error('Input argument no. 1 should be of type ''numeric'' and a matrix.')
                else
                    datapointer = libpointer('int32Ptr',varargin{1});
                    this.objectHandle = calllib('libgismo','gsMatrixInt_create_rc',varargin{1},size(varargin{1},1),size(varargin{1},2),datapointer);
                end
            elseif(nargin==2)
                this.objectHandle = calllib('libgismo','gsMatrixInt_create_rc',varargin{:});
            elseif(nargin==3)
                this.objectHandle = calllib('libgismo','gsMatrixInt_create_rcd',varargin{:});
            end            
        end

        function delete(this)
            calllib('libgismo','gsMatrixInt_delete',this.objectHandle);
        end

        % This function returns the address of the C++ pointer
        function varargout = ptr(this, varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            varargout{1} = this.objectHandle;
        end

        function disp(this, varargin)
            calllib('libgismo','gsMatrixInt_print',this.objectHandle);
        end


        function varargout = rows(this, varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsMatrixInt_rows',this.objectHandle);
        end

        function varargout = cols(this, varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            disp(nargout)
            [varargout{1:nargout}] = calllib('libgismo','gsMatrixInt_cols',this.objectHandle);
        end

        function varargout = data(this, varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            [varargout{1:nargout}] = calllib('libgismo','gsMatrixInt_data',this.objectHandle);
        end

        function varargout = asMatrix(this, varargin)
            if (nargin~=1 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            datapointer = calllib('libgismo','gsMatrixInt_data',this.objectHandle);
            setdatatype(datapointer,'int32Ptr',this.rows,this.cols);
            [varargout{1:nargout}] = datapointer.Value;
        end

    end
end