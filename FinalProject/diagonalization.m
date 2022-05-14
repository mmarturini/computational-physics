function out = diagonalization(in)
%Diagonalization restituisce una struttura contenente autovettori e
%autovalori di una Hamiltoniana con generico potenziale, costruita a
%partire da laplaciano con p=inf e con p finito

%%default

default.N = 1024;
default.distance = 0.01;
default.V = @(x) 0.5*x.^2;
default.type_of_oscillator = 'Harmonic';
default.BC = 'PBC';
default.p = 1;


%% input handling

if nargin == 0
    out = default;
    return;
end
if isempty(in)
    in = struct([]);
end

% fill from 'in' all new values, keeping the same ordering of default
in = recsetup(in,default);

%% main
switch in.type_of_oscillator    %optimal lattice spacing
    case 'Harmonic'  
        a = ((pi^2/2)/(in.V(in.N/2)))^(1/4);
    case 'Anharmonic'
        a = ((pi^2/2)/(in.V(in.N/2)))^(1/6);
    otherwise
        a = in.distance;
end

in.distance = a;
[H,x,Hmulti] = myhamiltonian1D(in.N,a,in.V,in.p,in.BC);

[u,e] = eig(H);
[u_1,e_1] = eig(Hmulti);
out.H.eigenvectors = u;
out.H.eigenvalues = diag(e);
out.Hmulti.eigenvectors = u_1;
out.Hmulti.eigenvalues = diag(e_1);
out.lattice = x;
out.in = in;
end
%% other functions

function s = recsetup(a,b)
    s = b;
    for fname = fieldnames(a)'
        if isstruct(a.(fname{1}))
            s.(fname{1}) = recsetup(a.(fname{1}),b.(fname{1}));
        else
            s.(fname{1}) = a.(fname{1});
        end
    end
end

