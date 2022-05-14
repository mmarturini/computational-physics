function out = tempEvol(in)
% temporal evolution of a 1-dimensional wave function with arbitrary potential.
% in = tempEvol returns the default setup as a struct.
% out = tempEvol(in) returns the structure out with field 'x','pdf' and field
% 'in' containing the corresponding setup, in addition to the final psi and time,
% in order to make it possible to continue the ongoing evolution.
% out = tempEvol([]) runs with default settings;

narginchk(0,1)

%% set defaults

default.number_of_points = 2e3;
default.reticular_step = 0.1;   
default.initial_time = 0;
default.time_step = 0.02;
default.time_span = 10;
default.number_of_steps = 5e3;
default.V = 0; 
default.BC = 'DBC';
default.psi0_is_gaussian = 1;
default.psi0params = [0 1 0]; %[-30 10 2]; % (x0, sigma, k0)
default.psi0 = [];
default.same_psi = 0;
default.scale_factor = 50;
default.operator_splitting = 0; % not recommended for big time_step
default.show_evolution = 1;

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

%% prep

N = in.number_of_points;
a = in.reticular_step;
x = -(N-1)/2:(N-1)/2; 
x = a*x';

pars = in.psi0params;
if ~in.same_psi || isempty(in.psi0)
    if in.psi0_is_gaussian
        psi = exp(1i*pars(3)*x).*exp(-(x-pars(1)).^2/pars(2)^2/2);
    else
        prompt = 'Enter psi in numeric vector form:\n ';
        psi = input(prompt);
        if length(psi) ~= N 
            error('numeric psi has wrong size')
        end
    in.initial_time = 0;    
    end
else
    psi = in.psi0;
end
psi = psi/norm(psi);
s = in.scale_factor;

t0 = in.initial_time;
dt = in.time_step;
ts = in.time_span;
if dt*in.number_of_steps < ts
    nt = in.number_of_steps;   
else
    nt = round(ts/dt);
end

if isnumeric(in.V)
    if isscalar(in.V)
        in.V = in.V*ones(size(x));
    end
    if length(in.V) ~= N
        error('numeric potential has wrong size')
    end
    V = in.V;
elseif isa(in.V,'function_handle')
    V = in.V(x);
else
    error('potential has wrong type')
end

%% main

f = figure; clf
set(f,'numbertitle','off','name','Temporal evolution')
plot(x,V);
hold on
hpsi = plot(x,s*abs(psi).^2);
ht = title(sprintf('t = %-5.2f, \t\t waiting for calculation...',t0));
drawnow;

if in.show_evolution
    if ~in.operator_splitting
        H = myhamiltonian1D(N,a,V,[],in.BC);
        out.Emean = real(psi'*H*psi);
        U = expm(-1i*H*dt);
        for j=1:nt
            psi = U*psi; 
            hpsi.YData = s*abs(psi).^2;
            ht.String = sprintf('t = %-5.2f',t0+j*dt);        
            drawnow; 
        end    
    else
        switch in.BC
            case 'DBC'    
                k = (pi/(N+1))*(1:N)';
            case 'NBC'
                k = (pi/N)*(0:N-1)';
            case 'PBC'
                n = floor(N/2);
                nn = floor((N-1)/2);
                k = (2*pi/N)*(-n:nn)';
            otherwise
                error('unknown type of boundary conditions)');
        end
        kinE = (k/a).^2/2;
        for j = 1:nt
            psi = exp(-1i*V*dt/2).*psi;
            switch in.BC
                case 'DBC'    
                    psi = idst(exp(-1i*kinE*dt).*dst(psi));
                case 'NBC'
                    psi = idct(exp(-1i*kinE*dt).*dct(psi));
                case 'PBC'
                    psi = isfft(exp(-1i*kinE*dt).*sfft(psi));
            end
            psi = exp(-1i*V*dt/2).*psi;
            hpsi.YData = s*abs(psi).^2;
            ht.String = sprintf('t = %-5.2f',t0+j*dt);
            drawnow; 
        end             
    end
else
    H = myhamiltonian1D(N,a,V,[],in.BC);
    out.Emean = real(psi'*H*psi);
    Utot = expm(-1i*H*nt*dt);
    psi = Utot*psi;
    hpsi.YData = s*abs(psi).^2;
    ht.String = sprintf('t = %-5.2f',t0+nt*dt);
    drawnow; 
end

in.psi0 = psi;
in.initial_time = t0+nt*dt;
out.in = in;
out.x = x;
out.pdf = abs(psi).^2;

% to obtain probabilities, plot(out.x,out.pdf) and then,for example
% sum(out.pdf(out.x > 2))

end

%% aux functions

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

