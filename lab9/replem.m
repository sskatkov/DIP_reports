function retval = repelem(element, varargin)

  if (nargin <= 1)
    error('repelem: Not enough input arguments')
    print_usage ();
    
  elseif (nargin == 2)
  
    v = varargin{1};

    if (isscalar (v))

      if (iscolumn (element))
        retval = element (:, ones (v, 1))'(:);

      elseif (isrow (element))
        retval = element (ones (v, 1), :)(:)';

      else
        error (['repelem: %gD Array objects require %g or more input ' ...
               'arguments, only %g given'], ...
                ndims (element), ndims (element) + 1, nargin);
        print_usage ();

      endif

    elseif (isvector (element)) && (isvector (v))

      if (numel (v) == numel (element))
        idx2 = prepareIdx (v);
        retval  = element (idx2);

      else # catch unequal element counts
        error (["repelem: varargin must either be scalar or have the same" ...
                " number of elements as the vector to be replicated"]);
        print_usage ();

      endif
      
    else # catch any arrays passed to element or varargin with nargin==2
      error (["repelem: when called with only two inputs they must be" ...
              " either scalars or vectors."]);
      print_usage ();

    endif
  
  
  elseif (nargin == 3)  #can simplify for known dimension count
    
    # avoid repeated function calls
    elsize = size (element);
    # 'numel' or 'length' faster than isvector in cellfun
    scalarv = (cellfun ('numel', varargin) == 1);
    nonscalarv = ~scalarv;
    
    ##INPUT CHECK

    #1:check that all varargin are either scalars or vectors, no arrays.
    # isvector gives true for scalars.
    # (Faster here with only two to avoid cellfun)
    if (~ (isvector (varargin{1}) && (isvector (varargin{2}))))
      error ("repelem: varargin must be all scalars or vectors");
      print_usage ();

    #2: check that the ones that are vectors have the right length.
    elseif (any (~ (cellfun ('numel', varargin (nonscalarv)) == elsize (nonscalarv))))
      error (["repelem: varargin(n) must either be scalar or have the same" ...
              " number of elements as the size of dimension n of the array" ...
              " to be replicated"]);
      print_usage ();

    endif 
    
    #Create index arrays to pass to element 
    ##(no slower passing to prepareIdx than checking and doing scalars directly)
    idx1 = prepareIdx (varargin{1}, elsize (1));
    idx2 = prepareIdx (varargin{2}, elsize (2));
  
    #the : at the end takes care of size(element)>2
    retval = element (idx1, idx2, :);

  else  #if (nargin > 3) **no need for elseif
  
    # avoid repeated function calls
    elsize = size (element);
    eldims = numel (elsize);
    vasize = nargin - 1; #numel(varargin);
    %maxDim = max(eldims,vasize);
    dimsWithBoth = min (eldims, vasize);
    nonscalarv   = ~cellfun (@isscalar, varargin);

    ## INPUT CHECK

    # 1: that they are all scalars or vectors. isvector gives true for scalars.
    if (~all (cellfun (@isvector, varargin)))
      error ("repelem: varargin must be all be scalars or vectors");
      print_usage ();

    # 2: catch any vectors thrown at trailing singletons, which should only have scalars
    elseif (max (find (nonscalarv)) > eldims)
      error ("repelem: varargin(n) for trailing singleton dimensions must be scalar");
      print_usage ();

    # 3: that the ones that are vectors have the right length.
    elseif (any (~ (cellfun ('numel', varargin (nonscalarv)) == elsize (nonscalarv))))
      error (["repelem: varargin(n) must either be scalar or have the same" ...
              " number of elements as the size of dimension n of the array" ...
              " to be replicated"]);
      print_usage ();

    endif 
    
    # first, preallocate idx which will contain index array to be put into element
    idx = cell (1, vasize);

    # use prepareIdx() to fill indices for each dimension that could be a scalar or vector
    idx (1:dimsWithBoth) = cellfun (@prepareIdx, varargin (1:dimsWithBoth), ...
                           num2cell (elsize (1:dimsWithBoth)), ...
                           'UniformOutput', false);

    # if there are more varargin inputs than element dimensions, input tests have 
    # verified they are just scalars, so add [1 1 1 1 1... 1] to those dims to 
    # perform concatenation in those dims.
    
    ###can cellfun check speed against for loop.
    if (vasize > eldims)
      idx (eldims + 1:vasize) = cellfun ('ones', {1}, ...
                                {varargin(eldims + 1:end){:}}, ...
                                'UniformOutput', false);
    endif
    
    # use completed idx to specify repitition of element values in all dimensions
    # trailing : will take care of any case where eldims > vasize.
    retval = element (idx{:}, :);
  
  endif

endfunction


function idx = prepareIdx (v, elsize_n)
# returns a row vector of indices prepared for replicating.

  if (isscalar (v))
    # will always return row vector
    idx = [1:elsize_n](ones (v, 1), :)(:)';

  else
  
    # works for row or column vector. idx2 output will be a row vector.
    # gets ending position for each element item
    idx_temp = cumsum (v);
    # row vector with enough space for output
    idx (1:idx_temp (end)) = 0;
    # sets starting position of each element to 1
    idx (idx_temp (1:end - 1) + 1) = 1;
    # sets starting position of each element to 1
    idx(1) = 1;
    # with prepared index
    idx = (find (v ~= 0))(cumsum (idx));

  endif
  
endfunction