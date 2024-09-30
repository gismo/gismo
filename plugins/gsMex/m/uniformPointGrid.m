function [pts] = uniformPointGrid(lower,upper,numPts)
%uniformPointGrid - Approximately uniformly spaced grid in every direction, with approximately numPts total points 
%
%Usage:
%  [pts] = uniformPointGrid(lower,upper,numPts)


%To avoid duplicate codes we call the underlying C++ method rather than implement it directly in MATLAB

  DIM = 2;
  if (~isa(lower,'numeric')  || numel(lower)~=DIM  || ...
      ~isa(upper,'numeric')  || numel(upper)~=DIM  || ...
      ~isa(numPts,'numeric') || numel(numPts)~=1   )
    error('Invalid input.')
  end
  pts = mex_uniformPointGrid(lower,upper,numPts);

end
