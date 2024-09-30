function varargout = gismo_map(pts, in, der)
  if ( iscell(pts) ) % create cartesian product points from cell infomation
      ndim = length(pts);
      npts = prod(cellfun(@length,pts));
      pts = cartesian_product_from_cell(pts);
  else
      [ndim, npts] = size(pts);
  end
  
  if (der == 0 || nargout > 1)
    F = in.eval(pts); % dimension rdim x npts
    varargout{1} = F;
  end
  if (der == 1 || nargout > 2)
    jac = reshape(in.deriv(pts),[],ndim,npts); % dimension rdim x ndim x npts
    if nargout == 1
        varargout{1} = jac;
    else % nargout == 2 or 3
        varargout{2} = jac;
    end
  end
  if (der == 2)
      rdim = length(in.eval(zeros(ndim,1)));
      hess = zeros(rdim,ndim,ndim,npts); % dimension rdim x ndim x ndim x npts
%       for dir = 1:rdim
%         hess(dir,:,:,:) = reshape(in.hess(pts,dir),ndim,ndim,npts);
%       end
    if nargout == 1
        varargout{1} = hess;
    elseif nargout == 3
        varargout{3} = hess;
    end
  end
end

function pts_aux = cartesian_product_from_cell(pts)
  % create cartesian product points from cell infomation
  s = cellfun(@length,pts);
  s_cell = cell(length(s),1);
  for ii = 1:length(s)
      s_cell{ii} = 1:s(ii);
  end
  x = cell(1,numel(s_cell));
  [x{:}] = ndgrid(s_cell{:});
  pts_aux = [];
  for ii=1:length(s)
      pts_aux = [pts_aux; reshape(pts{ii}(x{ii}),1,[])];
  end
end