function [c,lags] = xcorr(x,varargin)
%XCORR Cross-correlation function estimates.
%   C = XCORR(A,B), where A and B are length M vectors (M>1), returns
%   the length 2*M-1 cross-correlation sequence C. If A and B are of
%   different length, the shortest one is zero-padded. C will be a
%   row vector if A is a row vector, and a column vector if A is a
%   column vector.
%
%   XCORR produces an estimate of the correlation between two random
%   (jointly stationary) sequences:
%          C(m) = E[A(n+m)*conj(B(n))] = E[A(n)*conj(B(n-m))]
%   It is also the deterministic correlation between two deterministic
%   signals.
%
%   C = XCORR(A), where A is a length M vector, returns the length 2*M-1
%   auto-correlation sequence C. The zeroth lag of the output correlation
%   is in the middle of the sequence, at element M.
%
%   C = XCORR(A), where A is an M-by-N matrix (M>1), returns a large matrix
%   with 2*M-1 rows and N^2 columns containing the cross-correlation
%   sequences for all combinations of the columns of A; the first N columns
%   of C contain the delays and cross correlations using the first column
%   of A as the reference, the next N columns of C contain the delays and
%   cross correlations using the second column of A as the reference, and
%   so on.
%
%   C = XCORR(...,MAXLAG) computes the (auto/cross) correlation over the
%   range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags.
%   If missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCORR(...)  returns a vector of lag indices (LAGS).
%
%   XCORR(...,SCALEOPT), normalizes the correlation according to SCALEOPT:
%     'biased' - scales the raw cross-correlation by 1/M.
%     'unbiased' - scales the raw correlation by 1/(M-abs(lags)).
%     'normalized' or 'coeff' - normalizes the sequence so that the 
%                               auto-correlations at zero lag are 
%                               identically 1.0.
%     'none' - no scaling (this is the default).
%
%   % Example:
%   % Compute and plot the cross-correlation of two 16-sample
%   % exponential sequences
%
%   N = 16;
%   n = 0:N-1;
%   a = 0.84;
%   b = 0.92;
%   xa = a.^n;
%   xb = b.^n;
%   [r,lags] = xcorr(xa,xb);
%   stem(lags,r)
%
%   See also XCOV, CORRCOEF, CONV, COV.

%   Copyright 1988-2024 The MathWorks, Inc.

%   References:
%     S.J. Orfanidis, "Optimum Signal Processing. An Introduction"
%     2nd Ed. Macmillan, 1988.

narginchk(1,4);
if isnumeric(x)
    if isa(x,'uint64') || isa(x,'int64')
        error(message('MATLAB:xcorr:InvalidInputType'));
    end
elseif ~islogical(x)
    error(message('MATLAB:xcorr:InvalidInputType'));
end

[ySupplied,maxlagInput,scale] = ...
    matlab.internal.math.parseXcorrOptions(varargin{:});

% Transform the input so that computations can be performed on columns.
if size(x,1) == 1 && ~isscalar(x)
    % Shift the leading non-singleton dimension to the fore and make a
    % recursive call. Row vector input becomes a column vector, N-D
    % inputs have leading ones shifted out.
    [x,nshift] = shiftdim(x);
    if nargout == 2
        [c,lags] = xcorr(x,varargin{:});
    else
        c = xcorr(x,varargin{:});
    end
    c = shiftdim(c,-nshift);
    return
end

% Calculate the cross-correlation.
if ySupplied
    % Cross-correlation of two column vectors.
    if ~iscolumn(x)
        error(message('MATLAB:xcorr:MismatchedAB'));
    end
    y = varargin{1}(:); % y was validated to be a vector, make it a column.
    
    maxlagDefault = max(size(x,1),size(y,1)) - 1;
    if isempty(maxlagInput)
        maxlag = maxlagDefault;
    else
        maxlag = maxlagInput;
    end
    
    if isempty(x) || isempty(y)
        if isa(x,'single') || isa(y,'single')
            c = zeros(2*maxlag+1,1,'single');
        else
            c = zeros(2*maxlag+1,1);
        end
        lags = -maxlag:maxlag;
        return;
    else
        [c,lags] = xcorr(vA-mean(vA), vB-mean(vB), 'none'); % no scaling
        % scale to correlation coefficients yourself
        c = c ./ (std(vA) * std(vB));
    end
else
    % Perform all auto- and cross-correlations of the columns of x.
    narginchk(1,3);
    
    maxlagDefault = size(x,1) - 1;
    if isempty(maxlagInput)
        maxlag = maxlagDefault;
    else
        maxlag = maxlagInput;
    end
    
    if isempty(x)
        [~,n] = size(x);
        if isa(x,'single')
            c = zeros(2*maxlag+1,n^2,'single');
        else
            c = zeros(2*maxlag+1,n^2);
        end
        lags = -maxlag:maxlag;
        return;
    else
        % Peform the auto- and cross-correlations.
        c = autocorr(x,maxlag);
        
        % Scale the output.
        c = scaleOutput(scale,c,x);
    end
end

% Pad the output with zeros.
if maxlag > maxlagDefault
    zeropad = zeros(maxlag - maxlagDefault,size(c,2),"like",c);
    c = [zeropad; c; zeropad];
end

if nargout == 2
    lags = -maxlag:maxlag;
end

%--------------------------------------------------------------------------
function c = autocorr(x,maxlag)
% Compute all possible auto- and cross-correlations of the columns of a
% matrix input x. Output is clipped based on maxlag but not padded when
% maxlag >= size(x,1).
[m,n] = size(x);
mxl = min(maxlag,m - 1);

if n == 1
    % Autocorrelation of a column vector.
    if ~allfinite(x)
        c1 = conv(x,conj(flip(x)));
        c = c1((m-mxl):(m+mxl),1);
    else
        m2 = findTransformLength(m);
        X = fft(x,m2,1);
        Cr = abs(X).^2;
        if isreal(x)
            c1 = ifft(Cr,[],1,'symmetric');
        else
            c1 = ifft(Cr,[],1);
        end
        % Keep only the lags we want and move negative lags before positive
        % lags.
        c = [c1(m2 - mxl + (1:mxl),1); c1(1:mxl+1,1)];
    end
else
    % Auto- and cross-correlation of the columns of a matrix.
    if ~allfinite(x)
        ctr = 1;
        [m,n] = size(x);
        c1 = zeros(2*m-1,n);
        for k = 1:n
            for j = 1:n
                c1(:,ctr) = conv(x(:,k),conj(flip(x(:,j))));
                ctr = ctr + 1;
            end
        end
        c = c1((m-mxl):(m+mxl),:);
    else
        m2 = findTransformLength(m);
        X = fft(x,m2,1);
        % The number of rows of X is M2 (a power of 2), and X has N
        % columns. We want to perform all possible element-wise
        % multiplications of the columns of X by the columns of conj(X). To
        % do this efficiently, we use implicit expansion, asking it to
        % expand one operand in the second dimension and the other in the
        % third. Note that X should already be m2-by-n, but if we reshape
        % anyway via X(:,:), the code will be robust when operating on N-D
        % inputs.
        C = reshape(X,m2,1,n).*conj(X(:,:));
        % Call IFFT and force real output if x is real.
        if isreal(x)
            c1 = ifft(C,[],1,'symmetric');
        else
            c1 = ifft(C,[],1);
        end
        % c1 is M2-by-N-by-N.
        % Keep only the lags we want, and move the negative lags before the
        % positive lags. Also flatten the result to 2-D.
        c = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];
    end
end

%--------------------------------------------------------------------------
function c = crosscorr(x,y,maxlag)
% Compute cross-correlation for vector inputs. Output is clipped based on
% maxlag but not padded if maxlag >= max(size(x,1),size(y,1)).
nx = numel(x);
ny = numel(y);
m = max(nx,ny);
maxlagDefault = m-1;
mxl = min(maxlag,maxlagDefault);

if ~allfinite(x) || ~allfinite(y)
    c1 = conv(x,conj(flip(y)));
    if mxl <= maxlagDefault
        % Clip if maxlag is small. This would be a no-op for larger maxlag.
        % Account for 0 lag corresponding to c1(ny).
        c1 = c1(max(1,ny-mxl):(ny+min(mxl,nx-1)));
    end
    % Pad the head or tail if nx >= ny or nx < ny, respectively.
    % Note that we may need to both clip and pad: xcorr(1:5,1:3,3).
    c = zeros(2*mxl+1,1,"like",c1);
    offset = max(0,mxl-ny+1);
    c(offset+(1:numel(c1))) = c1;
else
    m2 = findTransformLength(m);
    X = fft(x,m2,1);
    Y = fft(y,m2,1);
    if isreal(x) && isreal(y)
        c1 = ifft(X.*conj(Y),[],1,'symmetric');
    else
        c1 = ifft(X.*conj(Y),[],1);
    end
    % Keep only the lags we want and move negative lags before positive
    % lags.
    c = [c1(m2 - mxl + (1:mxl)); c1(1:mxl+1)];
end

%--------------------------------------------------------------------------
function m = findTransformLength(m)
m = 2*m;
while true
    r = m;
    for p = [2 3 5 7]
        while (r > 1) && (mod(r, p) == 0)
            r = r / p;
        end
    end
    if r == 1
        break;
    end
    m = m + 1;
end

%--------------------------------------------------------------------------
function c = scaleOutput(scale,c,x,y)
% Scale correlation as specified.
if strcmp(scale,'none')
    return
end

ySupplied = nargin == 4;
m = size(x,1);
if ySupplied && (m ~= size(y,1))
    error(message('MATLAB:xcorr:NoScale'));
end

if strcmp(scale,'biased')
    % Scales the raw cross-correlation by 1/M.
    c = c./m;
elseif strcmp(scale,'unbiased')
    % Scales the raw correlation by 1/(M-abs(lags)).
    L = (size(c,1) - 1)/2;
    scaleUnbiased = (m - abs(-L:L)).';
    scaleUnbiased(scaleUnbiased <= 0) = 1;
    c = c./scaleUnbiased;
else % 'normalized'/'coeff'
    % Normalizes the sequence so that the auto-correlations
    % at zero lag are identically 1.0.
    if ySupplied
        % Compute autocorrelations at zero lag.
        % scale = norm(x)*norm(y) is numerically superior but slower.
        cxx0 = sum(abs(x).^2);
        cyy0 = sum(abs(y).^2);
        scaleCoeffCross = sqrt(cxx0*cyy0);
        c = c./scaleCoeffCross;
    elseif size(c,2) == 1
        % Autocorrelation of a vector. Normalize by c[0].
        mid = (size(c,1) + 1)/2; % row corresponding to zero lag.
        c = c./c(mid);
    else
        % Compute the indices corresponding to the columns that are
        % autocorrelations.
        [~,n] = size(x);
        % Note that size(c,2) = n^2.
        kvec = 1:n+1:n*n; % a 1-by-n row vector
        % kvec is an index vector such that for an n-by-n matrix A,
        % A(kvec(j)) = A(j,j).
        mid = (size(c,1) + 1)/2; % row index corresponding to zero lag
        trow = sqrt(c(mid,kvec)); % a 1-by-n row vector
        tmat = trow.'*trow; % an n-by-n matrix, tmat(i,j) = trow(i)*trow(j)
        scaleCoeffAuto = tmat(:).'; % a 1-by-n^2 row vector
        % The autocorrelations at zero-lag are normalized to one.
        c = c./scaleCoeffAuto;
    end
end
