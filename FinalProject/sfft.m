function [tf,k] = sfft(f,x)
% symmetric fast-fourier  tf = sfft(f);
% acts columnwise unless f is a row vector;
% normalization as in standard fft, that is tf contains N times the 
% coefficients of the Fourier series where the the original function 
% has the same dimensions of its Fourier coefficients./
%
% [tf,k] = sfft(f,x) returs also the properly normalized wavenumbers,
% using only the lattice spacing a=x(2)-x(1) of uniform array x 
% representing the discretized real segment. In this case tf is normalized
% as in the full space Fourier transform, i.e. [tf]=[x]*[f], with no
% factor in the f --> tf map and the (1/2/pi) the tf --> f map.

narginchk(1,2)

[N,nc] = size(f);
flag = isrow(f);
if flag
    f = f.'; N = nc; nc = 1 ;
end

p = floor((N-1)/2);
n = [0:p,p+1-N:-1]';
tf = exp(1i*pi*(1-1/N)*n)*ones(1,nc) .* fft(f);
tf = tf([p+2:N,1:p+1],:);

if flag
    tf = tf.';
end

k = [];
if nargin == 2 && nargout == 2
    a = x(2)-x(1);
    pp = floor(N/2);
    k = (2*pi/N/a)*(-pp:p);
    if iscolumn(x)
        k = k.';
    end
    tf = a*tf;
end
    
    

