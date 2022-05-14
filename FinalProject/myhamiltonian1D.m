function [H,x,Hmulti] = myhamiltonian1D(N,a,V,p,BC) 
%myhamiltonian1D riceve come input: 
%N = numero di punti del reticolo
%a = passo reticolare 
%V = potenziale
%2p+2 = ordine di sviluppo del laplaciano
%BC = condizioni al contorno
%Restituisce come output:
%H = hamiltoniana costruita da laplciano pieno
%Hmulti = hamiltoniana costruita da laplaciano multidiagonale
%x = reticolo

% the grid 
x = -(N-1)/2:(N-1)/2;
x = a*x';

if isnumeric(V)
    if isscalar(V)
        V = V*ones(size(x));
    end
    if length(V) ~= N
        error('numeric potential has wrong size')
    end
elseif isa(V,'function_handle')
    V = V(x);
else
    error('potential has wrong type')
end

switch BC
    case 'DBC'
       L = mylaplacian1D(N,inf,'DSTI');
    case 'NBC'
       L = mylaplacian1D(N,inf,'DCTII');
    case 'PBC'
       L = mylaplacian1D(N,inf,'PBC');
    otherwise
        error('unknown type of boundary conditions')
end

H = L/2/a^2 + diag(V);
if ~isempty(p) 
    Hmulti =  full(mylaplacian1D(N,p,BC))/2/a^2 + diag(V); 
end

end