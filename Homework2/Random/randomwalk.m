function out = randomwalk(in)
% n-dimensional random walk
% in = randomwalk returns the default setup as a struct.
% out = randomwalk(in) returns the structure out with field 
% 'final_positions' and field 'in' containing the corresponding setup.
% out = randomwalk([]) runs with default settings;


% TODO implement variable step update algorithm?

narginchk(0,1)

%% set defaults (2-dim case)

default.number_of_particles = 1e4;
default.boundary_conditions = [0 0 0 0];
default.initial_positions = @(N)0.01*randn(N,2);
default.drift_field = @(t,x)zeros(size(x));
default.random_jumps = @(t,x)randn(size(x));
default.time_step = 1e-4; % for fixed step algorithms
default.step_tolerance = 1e-6; % for variable step algorithms
default.time_span = 1;
default.step_algorithm = 'euler';
default.space_is_phase_space = 0;
default.plot.viewbox = [-1 1 -1 1];
default.plot.figure_number = 5;
default.plot.figure_position = [1400 600 560 420];
default.plot.frame_skips = 0;
default.plot.uicontrols = true;
default.hist.varpair = 'xr';  % '' 'xy' 'xz' 'xr' ....
default.hist.nbin = 30;
default.hist.maxheight = [5 5];
default.parallel_procs = 1;
default.save_walks = false; % careful with memory usage

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

% fill short-named variables and perform some consistency check
vbox = in.plot.viewbox;
dim = length(vbox)/2;
bc = zeros(1,2*dim); % default free boundary conditions
bc(1:length(in.boundary_conditions)) = in.boundary_conditions;

N = in.number_of_particles;
nproc = max(1,in.parallel_procs);
Nperproc = floor(N/nproc);
if N ~= nproc*Nperproc
    warning('number of particles not divisible by number of processes')
    N = nproc*Nperproc;
end
if isnumeric(in.initial_positions)
    x = in.initial_positions;
else
    x = in.initial_positions(N);
end
if nproc > 1
    x = reshape(x,Nperproc,nproc,dim);
    x = permute(x,[1 3 2]);
end

B = in.drift_field;
rj = in.random_jumps;
dt = in.time_step;
nt = round(in.time_span/dt);
if dim ~= 2 && dim ~= 3
    skip = inf; % do not plot in dimensions other that 2 or 3
else
    skip = in.plot.frame_skips;
end

update = @update0;
if in.space_is_phase_space
    if mod(dim,2) ~= 0
        error('a classical phase space must have even dimensionality')
    end
    update = @update1;
end

onestep = in.step_algorithm;

%% main

% handle graphics
run.Value = 1;
stop.Value = 0;
if skip < nt && nproc == 1
    hf = figure(in.plot.figure_number);
    clf(in.plot.figure_number,'reset');
    set(hf,'numbertitle','off','name','Random walk')
    set(hf,'position',in.plot.figure_position);
    ha = axes('outerposition',[0 .02 1 1]);
    if dim == 2
        h = plot(ha,x(:,1),x(:,2),'.');
    else
        h = plot3(ha,x(:,1),x(:,2),x(:,3),'.');
    end
    box('on'); axis(vbox); axis square; shg
    if ~isempty(in.hist.varpair)
        nbin = in.hist.nbin;
        set(hf,'position',in.plot.figure_position + [0 -250 0 250]);
        ha.OuterPosition = [0.1 .3 0.8 0.8];
        hb = axes('outerposition',[0 0.05 0.55 0.375]);
        hc = axes('outerposition',[0.45 0.05 0.55 0.38]);
        [xb,binb] = histplot(in.hist.varpair(1),x);
        [xc,binc] = histplot(in.hist.varpair(2),x);
        hxb = histogram(hb,xb,binb,'Normalization','pdf');
        hxc = histogram(hc,xc,binc,'Normalization','pdf');
        hc.YTickLabel = [];
        hb.YLim = [0 in.hist.maxheight(1)];
        hc.YLim = [0 in.hist.maxheight(2)];
    end
    ttl = title(sprintf('t = %-8.1f',0));
    if in.plot.uicontrols 
        run = uicontrol('style','toggle','string','start','value',0,...
            'units','normalized','position',[.375 .01 .1 .05]);
        stop = uicontrol('style','toggle','string','stop','value',0,...
            'units','normalized','position',[.55 .01 .1 .05]);
    end
    drawnow
end

if nproc > 1 && isfinite(nt)
    parfor jproc = 1:nproc
        for jt = 1:nt
            t = jt*dt;
            x(:,:,jproc) = update(t,x(:,:,jproc)); %#ok<PFBNS>
        end
    end
else
    jt = 0;
    while ~stop.Value && jt < nt
        if ~run.Value
            pause(0.2)
            if strcmp(run.String,'pause')
                set(run,'string','resume','callback','run.Value=1;')
            end
        else
            jt = jt+1;
            t = jt*dt;
            x = update(t,x);
            if skip < nt
                if mod(jt+1,skip+1) == 0
                    h.XData = x(:,1);
                    h.YData = x(:,2);
                    if dim == 3
                        h.ZData = x(:,3);
                    end
                    if ~isempty(in.hist.varpair)
                        hxb.Data = histplot(in.hist.varpair(1),x);
                        hxc.Data = histplot(in.hist.varpair(2),x);
                    end
                    ttl.String = sprintf('t = %-5.2f',t);
                    drawnow
                end
                if in.plot.uicontrols
                    set(run,'string','pause','callback','run.Value=0;')
                end
            end
        end
    end
end

if nproc > 1
    x = permute(x,[1 3 2]);
    x = reshape(x,N,dim);
end

out.final_positions = x;
out.in = in;

if skip < nt && nproc == 1 && in.plot.uicontrols 
    set(stop,'string','close','value',0,'callback','close(gcf)')
end

%% nested funcs

    function [u,bin] = histplot(w,x)
        if w == 'x'
            u = x(:,1);
            bin = linspace(1.25*vbox(1),1.25*vbox(2),nbin);
        elseif w == 'y'
            u = x(:,2);
            bin = linspace(1.25*vbox(3),1.25*vbox(4),nbin);
        elseif w == 'z' && dim == 3
            u = x(:,1);
            bin = linspace(1.25*vbox(5),1.25*vbox(6),nbin);
        elseif w == 'r'
            u = sqrt(sum(x.^2,2));
            dbox = sqrt(sum(vbox(2:2:end).^2));
            bin = linspace(0,1.25*dbox,nbin);
        end
    end

    function x = update0(t,x)
        xnew = feval(onestep,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x);
        %x = xnew + sqrt(dt)*(rj(t,x) + rj(t+dt,xnew))/2;
        for jdim = 1:dim
            jo = 2*jdim-1;
            je = 2*jdim;
            if bc(jo)
                I = x(:,jdim) < vbox(jo);
                x(I,jdim) = -x(I,jdim) + 2*vbox(jo);
            end
            if bc(je)
                I = x(:,jdim) > vbox(je);
                x(I,jdim) = -x(I,jdim) + 2*vbox(je);
            end
        end
    end

    function x = update1(t,x)
        xnew = feval(onestep,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x);
        %x = xnew + sqrt(dt)*(rj(t,x) + rj(t+dt,xnew))/2;
        for jdim = 1:dim/2
            jo = 2*jdim-1;
            je = 2*jdim;
            if bc(jo)
                I = x(:,jdim) < vbox(jo);
                x(I,jdim) = -x(I,jdim) + 2*vbox(jo);
                x(I,jdim+dim/2) = -x(I,jdim+dim/2);
            end
            if bc(je)
                I = x(:,jdim) > vbox(je);
                x(I,jdim) = -x(I,jdim) + 2*vbox(je);
                x(I,jdim+dim/2) = -x(I,jdim+dim/2);
            end
        end
    end



end

%% aux functions

function yp = euler(F,t,y,h) %#ok<DEFNU>
    yp = y + h*F(t,y);
end

function yp = rk4(F,t,y,h) %#ok<DEFNU>
    s1 = F(t,y);
    s2 = F(t+h/2,y+h/2*s1);
    s3 = F(t+h/2,y+h/2*s2);
    s4 = F(t+h,y+h*s3);
    yp = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
end

