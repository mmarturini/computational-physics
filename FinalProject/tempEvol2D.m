function out = tempEvol2D(in)
% temporal evolution of a 2-dimensional wave function with arbitrary potential.
% in = tempEvol returns the default setup as a struct.
% out = tempEvol(in) returns the structure out with field 'x','y','pdf' and field
% 'in' containing the corresponding setup, in addition to the final psi and time,
% in order to make it possible to continue the ongoing evolution.
% out = tempEvol([]) runs with default settings;

narginchk(0,1)

%% set defaults

default.number_of_points_x = 200;
default.number_of_points_y = 200;
default.reticular_step = 0.05;   
default.initial_time = 0;
default.time_step = 0.05;
default.time_span = 10;
default.number_of_steps = 5e3;
default.V = @(x,y)x.^2/2+y'.^2/2;
%@(x,y)50*double((x+1).^2<0.04)*double(abs(y.^2-0.04)>0.02)';
default.BC_x = 'DBC';
default.BC_y = 'DBC';
default.psi0_is_gaussian = 1;
default.psi0_gaussian_params = [-2.5 2.5 1 1 0 0];  %[x0, y0, sigmax, sigmay, k0x, k0y]
%[-3 0 0.6 0.6 30 0];
default.psi0 = [];
default.same_psi = 1;
default.scale_factor = 2e4;
default.hZ = [0 20];
default.plot_color = 'parula';
default.operator_splitting = 1; % highly recommended for small time_step and big Nx, Ny
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

Nx = in.number_of_points_x;
Ny = in.number_of_points_y;
a = in.reticular_step;
x = -(Nx-1)/2:(Nx-1)/2; 
y = -(Ny-1)/2:(Ny-1)/2; 
x = a*x';
y = a*y';
d_mat = [Nx,Ny];
d_vec = [Nx*Ny,1];

if ~in.same_psi || isempty(in.psi0)
    if in.psi0_is_gaussian
        pars = in.psi0_gaussian_params;    
        psix = exp(1i*pars(5)*x).*exp(-(x-pars(1)).^2/pars(3)^2/2);
        psiy = exp(1i*pars(6)*y).*exp(-(y-pars(2)).^2/pars(4)^2/2);
        psi = psix*psiy';       
    else
        prompt = 'Enter psi in numeric matrix form:\n ';
        psi = input(prompt);
        if size(psi,1) ~= Nx || size(psi,2) ~= Ny
            error('numeric psi has wrong size')
        end
    end
    in.initial_time = 0;
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
        in.V = in.V*ones(Nx,Ny);
    end
    if size(in.V,1) ~= Nx || size(in.V,2) ~= Ny
        error('numeric potential has wrong size')
    end
    V = in.V;
elseif isa(in.V,'function_handle')
    V = in.V(x,y);
else
    error('potential has wrong type')
end

%% main

f = figure(1); clf
set(f,'numbertitle','off','name','2D Temporal evolution')
mesh(y,x,V);
if ~isempty(in.hZ)
    zlim(in.hZ);
end
if ~isempty(in.plot_color)
    colormap(in.plot_color);
end
hold on
hpsi = surf(y,x,s*abs(psi).^2);
hpsi.EdgeColor = 'none';
ht = title(sprintf('t = %-5.3f, \t\t waiting for calculation...',t0));

if in.show_evolution
    if ~in.operator_splitting
        psi = reshape(psi,d_vec);
        H = hamiltonian2D(Nx,Ny,a,V,in.BC_x,in.BC_y);
        U = expm(-1i*H*dt);
        for j=1:nt
            psi = U*psi; 
            hpsi.ZData = s*abs(reshape(psi,d_mat)).^2;
            ht.String = sprintf('t = %-5.3f',t0+j*dt);        
            drawnow; 
        end    
        psi = reshape(psi,d_mat);
    else
        kx = find_k(Nx,in.BC_x);
        ky = find_k(Ny,in.BC_y);
        kinx_exp = exp(-1i*dt*(kx/a).^2/2);
        kiny_exp = exp(-1i*dt*(ky/a).^2/2);
        V_exp = exp(-1i*V*dt/2);
        for j=1:nt           
            psi = V_exp.*psi;
            switch in.BC_x
                case 'DBC'
                    psi = idst(kinx_exp.*dst(psi));
                case 'NBC'
                    psi = idct(kinx_exp.*dct(psi));
                case 'PBC'
                    psi = isfft(kinx_exp.*sfft(psi));
            end
            switch in.BC_y
                case 'DBC'
                    psi = idst(kiny_exp.*dst(psi.')).';
                case 'NBC'
                    psi = idct(kiny_exp.*dct(psi.')).';
                case 'PBC'
                    psi = isfft(kiny_exp.*sfft(psi.')).';
            end
            psi = V_exp.*psi;
            hpsi.ZData = s*abs(psi).^2;
            ht.String = sprintf('t = %-5.3f',t0+j*dt);        
            drawnow; 
        end                 
    end
else   
    psi = reshape(psi,d_vec);
    H = hamiltonian2D(Nx,Ny,a,V,in.BC_x,in.BC_y);
    Utot = expm(-1i*H*nt*dt);
    psi = Utot*psi;
    psi = reshape(psi,d_mat);
    hpsi.ZData = s*abs(psi).^2;
    ht.String = sprintf('t = %-5.3f',t0+nt*dt);
    drawnow; 
end

in.psi0 = psi;
in.initial_time = t0+nt*dt;
out.in = in;
out.x = x;
out.y = y;
out.pdf = abs(psi).^2;
out.Vmean = Vmean;
out.rmean = rmean;
out.t = t0+dt:dt:nt*dt;

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

function L = laplacian1D(N,BC)
switch BC
    case 'DBC'
        k = (pi/(N+1))*(1:N)';
        L = idst(diag(k.^2)*dst(eye(N)));
    case 'NBC'
        k = (pi/N)*(0:N-1)';
        L = idct(diag(k.^2)*dct(eye(N)));
    case 'PBC'
        n = floor(N/2);
        nn = floor((N-1)/2);
        k = (2*pi/N)*(-n:nn)';
        L = isfft(diag(k.^2)*sfft(eye(N)));
    otherwise
        error('unknown type of boundary conditions')
end
L = real(L+L')/2;
end

function H = hamiltonian2D(Nx,Ny,a,V,BC_x,BC_y)
    Lx = laplacian1D(Nx,BC_x);
    Ly = laplacian1D(Ny,BC_y);
    L = kron(Ly,speye(Nx)) + kron(speye(Ny),Lx);
    H = L/2/a^2+diag(reshape(V,[Nx*Ny,1]));
end

function k = find_k(N,BC)
    switch BC
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
end


