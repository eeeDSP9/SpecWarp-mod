function fm_warp(input_file, output_file)    

% Load a speech waveform
    [d,sr] = audioread(input_file);
    [a,g,e] = lpcfit(d,20);
    alpha = -0.1;
    [bhat, ahat]  = warppoles(a, alpha);

    dw = filter(bhat(1,:), 1, lpcsynth(ahat, g, e));

    audiowrite(output_file,dw,sr)
end
function [a,g,e] = lpcfit(x,p,h,w,ov)
% [a,g,e] = lpcfit(x,p,h,w,ov)  Fit LPC to short-time segments
%    x is a stretch of signal.  Using w point (2*h) windows every 
%    h points (128), fit order p LPC models.  Return the successive 
%    all-pole coefficients as rows of a, the per-frame gains in g 
%    and the residual excitation in e.
%    ov nonzero selects overlap-add of window-length
%    residuals, otherwise successive hop-sized residuals are concatenated
%    for independent near-perfect reconstruction with lpcsynth.
%    (default is 1)
% 2001-02-25 dpwe@ee.columbia.edu $Header: /homes/dpwe/matlab/columbiafns/RCS/lpcfit.m,v 1.1 2004/03/30 20:55:52 dpwe Exp $

if nargin < 2
  p = 12;
end
if nargin < 3
  h = 128;
end
if nargin < 4
  w = 2*h;
end
if nargin < 5
  ov = 1;
end

if (size(x,2) == 1)
  x = x';  % Convert X from column to row
end

npts = length(x);

nhops = floor(npts/h);

% Pad x with zeros so that we can extract complete w-length windows
% from it
x = [zeros(1,(w-h)/2),x,zeros(1,(w-h/2))];

a = zeros(nhops, p+1);
g = zeros(nhops, 1);
if ov == 0
  e = zeros(1, npts);
else
  e = zeros(1, (nhops-1)*h+w);
end

% Pre-emphasis
pre = [1 -0.9];
x = filter(pre,1,x);

for hop = 1:nhops

  % Extract segment of signal
  xx = x((hop - 1)*h + [1:w]);
  % Apply hanning window
  wxx = xx .* hanning(w)';
  % Form autocorrelation (calculates *way* too many points)
  rxx = xcorr(wxx);
  % extract just the points we need (middle p+1 points)
  rxx = rxx(w+[0:p]);
  % Setup the normal equations
  R = toeplitz(rxx(1:p));
  % Solve for a (horribly inefficient to use full inv())
  an = inv(R)*rxx(2:(p+1))';
  % Calculate residual by filtering windowed xx
  aa = [1 -an'];
  if ov == 0
    rs = filter(aa, 1, xx((w-h)/2 + [1:h]));
  else
    rs = filter(aa,1,wxx);
  end
  G = sqrt(mean(rs.^2));
  % Save filter, gain and residual
  a(hop,:) = aa;
  g(hop) = G;
  if ov == 0
    e((hop - 1)*h + [1:h]) = rs'/G;
  else
    e((hop - 1)*h + [1:w]) =  e((hop - 1)*h + [1:w]) + rs/G;
  end
end

% Throw away first (win-hop)/2 pts if in overlap mode
% for proper synchronization of resynth
if ov ~= 0
  e = e((1+((w-h)/2)):end);
end

end
function d = lpcsynth(a,g,e,h,ov)
% d = lpcsynth(a,g,e,h,ov)   Resynthesize from LPC representation
%    Each row of a is an LPC fit to a h-point (non-overlapping) 
%    frame of data.  g gives the overall gains for each frame and 
%    e is an excitation signal (if e is empty, white noise is used; 
%    if e is a scalar, a pulse train is used with that period).
%    ov nonzero selects overlap-add of reconstructed 
%    windows, else e is assumed to consist of independent hop-sized 
%    segments that will line up correctly without cross-fading
%    (matching the ov option to lpcfit; default is ov = 1).
%    Return d as the resulting LPC resynthesis.
% 2001-02-25 dpwe@ee.columbia.edu $Header: /homes/dpwe/matlab/columbiafns/RCS/lpcsynth.m,v 1.1 2004/03/30 20:56:04 dpwe Exp $

if nargin < 3
  e = [];
end
if nargin < 4
  h = 128;
end
if nargin < 5
  ov = 1;
end

w = 2*h;

[nhops,p] = size(a);

npts = nhops*h;
% Excitation needs extra half-window at the end if in overlap mode
nepts = npts + ov*(w-h);

if length(e) == 0
  e = randn(1,nepts);
elseif length(e) == 1
  pd = e;
  e = sqrt(pd) * (rem(1:nepts,pd) == 0);
else
  nepts = length(e);
  npts = nepts - ov*(w-h);
end

% Try to make sure we don't run out of e (in ov mode)
e = [e, zeros(1, w)];

d = zeros(1,npts);

for hop = 1:nhops
  
  hbase = (hop-1)*h;
  
  oldbit = d(hbase + [1:h]);
  aa = a(hop,:);
  G = g(hop);
  if ov == 0 
    newbit = G*filter(1, aa, e(hbase + [1:h]));
  else
    newbit = G*filter(1, aa, e(hbase + [1:w]));
  end
  if ov == 0
    d(hbase + [1:h]) = newbit;
  else
    d(hbase + [1:w]) = [oldbit, zeros(1,(w-h))] + (hanning(w)'.*newbit); 
  end
  
end

% De-emphasis (must match pre-emphasis in lpcfit)
pre = [1 -0.9];
d = filter(1,pre,d);
end
function [B,A] = warppoles(a,alpha)
%  [B,A] = warppoles(a,alpha)  warp an all-pole polynomial by substitution
%    Warp all-pole polynomials defined by rows of a by the first-order 
%    warp factor alpha.  Negative alpha shifts poles up in frequency.
%    Output polynomials have zeros too, hence B and A.
% 2003-12-10 dpwe@ee.columbia.edu

% Construct z-hat^-1 polynomial
d = [-alpha 1];
c = [1 -alpha];

[nrows,order] = size(a);

A = zeros(nrows, order);
B = zeros(nrows, order);

B(:,1) = a(:,1);
A(:,1) = ones(nrows,1);

dd = d;
cc = c;

% This code originally mapped zeros.  I adapted it to map
% poles just by interchanging b and a, then swapping again at the 
% end.  Sorry that makes the variables confusing to read.

for n = 2:order

  for row = 1:nrows

    % add another factor to num, den
    B(row,1:order) = conv(B(row,1:(order-1)), c);

  end

  % accumulate terms from this factor
  B(:,1:length(dd)) = B(:,1:length(dd)) + a(:,n)*dd;
    
  dd = conv(dd, d);
  cc = conv(cc, c);

end
    
% Construct the uniform A polynomial (same for all rows)
AA = 1;
for n = 2:order
  AA = conv(AA,c);
end
A = repmat(AA, nrows, 1);
  

% Exchange zeros and poles
T = A; A = B; B = T;
end
