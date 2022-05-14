function [T,k] = Tcoeff(V,N,a,D2)
% [T,k] = Tcoeff(V,N,L) returns the transmission coefficient T 
% over a set of positive wavenumbers k for the symmetric potential V 
% using a grid with N point and lattice spacing and the Laplacian matrix
% D2, which defaults to laplacian1D(N,inf,'PBC').

narginchk(3,4)
if nargin < 4
    D2 = mylaplacian1D(N,inf,'PBC');
end

% the grid 
x = -(N-1)/2:(N-1)/2;
x = a*x';

% input check
if isnumeric(V) 
    if length(V) ~= N
        error('numeric potential has wrong size')
    end
elseif isa(V,'function_handle')
    V = V(x);
else
    error('potential has wrong type')
end
if ~isequal(V,flip(V))
    error('potential is not symmetric')
end
if size(D2,1) ~= N
        error('Laplacian matrix has wrong size')
end
% the hamiltonian
H = D2/2/a^2 + diag(V);

% diagonalization and phase shifts
E = sort(eig(H));
E(E<0) = []; % keep only scattering part
ke = sqrt(2*E(1:2:end));
ko = sqrt(2*E(2:2:end));
xie = pi*(0:length(ke)-1)' - N*a*ke/2;
xio = pi*(1:length(ko))' - N*a*ko/2;
        
% interpolation
k = linspace(ke(1),max(ke(end),ko(end)),1e4)';
fxi = griddedInterpolant(ke,xie,'pchip');
xie = fxi(k);
fxi = griddedInterpolant(ko,xio,'pchip');
xio = fxi(k);

% transmission coefficient
T = cos(xie-xio).^2;



