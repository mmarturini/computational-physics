function out = brownian1D(in)
% One-dimensional Brownian motion as two-dim random walk
% in = brownian1D returns default setting struct;
% out = brownian1D(in) returns the structure out with fields
% 'final_positions', 'final_velocities' and 'in', which contains the 
% setting struct used.
% out = brownian1D([]) runs with default settings;

% N.B.: this routine calls randomwalk.m internally

narginchk(0,1)

%% set defaults

default.number_of_particles = 5e4;
default.boundary_conditions = [1 0];
default.initial_positions = @(N)10*rand(N,1);
default.initial_velocities = @(N)zeros(N,1);
default.acceleration_field = @(t,x)-ones(size(x));
default.friction_coeff = 1;
default.kT_over_mass = 1;
default.random_kicks = @(t,x)randn(size(x));
default.time_step = 0.001; % for fixed step algorithms
default.time_span = 30;
default.step_algorithm = 'rk4';
default.plot.viewbox = [0 6 -5 5];
default.plot.frame_skips = 19;
default.hist.varpair = '';
default.parallel_procs = 1;

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

%% main: adapt input and call randomwalk

gamma = in.friction_coeff;
a = in.acceleration_field;
D_v = in.kT_over_mass*gamma;
rk = in.random_kicks;
x = in.initial_positions;
v = in.initial_velocities;

default = randomwalk;
% fill values from corresponding names
in2 = recsetup(in,default);

% adjust settings
in2.boundary_conditions(3:4) = [0 0];
if isnumeric(x)
    in2.initial_positions = [x v];
else
    in2.initial_positions = @(N)[x(N) v(N)];
end
in2.drift_field = @(t,x)[x(:,2), -gamma*x(:,2) + a(t,x(:,1))]; 
in2.random_jumps = @(t,x)[zeros(size(x,1),1), sqrt(2*D_v)*rk(t,x(:,2))];
in2.space_is_phase_space = 1;

o = randomwalk(in2);
out.final_positions = o.final_positions(:,1); 
out.final_velocities = o.final_positions(:,2);
out.in = in;

end
