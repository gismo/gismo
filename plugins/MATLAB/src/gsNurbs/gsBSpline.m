% @file gsBSpline.m
%
%    @brief Matlab wrapper for gsBSpline class
%
%    This file is part of the G+Smo library.
%
%    This Source Code Form is subject to the terms of the Mozilla Public
%    License, v. 2.0. If a copy of the MPL was not distributed with this
%    file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%    Author(s): H.M. Verhelst

classdef gsBSpline < gsGeometry

    methods(Access = public)

        % Constructor - Create a new C++ class instance 
        function this = gsBSpline(varargin)
            %gsBSpline - construct a gsBSpline object
            %
            %Usage:
            %  TODO
            %
            %Input:
            %  TODO
            %
            %Output:
            %  TODO

            if (~nargin==2 || nargout>1)
                error('Invalid number of input and/or output arguments.')
            end
            
            if ( isa(varargin{1},'gsBSplineBasis') && ismatrix(varargin{2}) )
                assert(varargin{1}.size()==size(varargin{2},1), 'Size mismatch between basis and coefficients.')
                coefs = EigenMatrix(size(varargin{2},1),size(varargin{2},2),varargin{2});
                this.objectHandle = calllib('libgismo','gsBSpline_create',varargin{1}.ptr(),coefs.ptr());
            else
                error('Input argument no. 1 should be of type ''gsBSplineBasis'' and argument no. 2 should be a matrix.')
            end
        end
    end
end
