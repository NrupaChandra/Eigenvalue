% plot_contour_deformation_v4_7.m
%
% plot_contour_deformation.m, with ONE thing changed: the JSON loader.
% Everything below the "Load JSON" heading is the fix; the driver and the
% plotting functions are unchanged from your original.
%
% WHY THE ORIGINAL FAILS on v4.6 / v4.7 logs
%   briggsv4.x.jl writes the iteration-1 entry from a different Dict than
%   every later entry: 15 keys for iteration 1, 24 for the rest in v4.7
%   (19 in v4.6).  jsondecode() builds a STRUCT ARRAY only when every
%   object in the JSON array has the same field names.  One odd entry and
%   it returns a 501x1 CELL ARRAY of structs instead.  jsonData(k) is then
%   a 1x1 cell, and
%       entry.iteration
%   throws "Dot indexing is not supported for variables of this type."
%   on the very first pass of the loop, so nothing is plotted at all.
%   v4.5 wrote both entries from the same 15-key Dict, which is why
%   contour_video_v4.5.avi exists and v4.6 has no video.
%
% The loader below normalises both shapes to a cell array and reads fields
% defensively, so it works on v4.5, v4.6, v4.7 and anything later that adds
% more diagnostics.

%% ---------------------------------------------------------------- settings
JSON_FILE  = 'contour_iteration_v4.7.json';
OUT_FILE   = 'contour_video_v4.7.mp4';
ITER_START = 1;      % first entry to render (1 = keep original behaviour)
ITER_STRIDE= 1;      % 1 = every entry.  4 or 8 renders a much shorter video.

data = read_contour_integration_json(JSON_FILE);

%% Video axis limits
numIterations = length(data);
omega_real = [];
omega_imag = [];
L_real = [];
L_imag = [];
F_real = [];
F_imag = [];
alphaU_real = [];
alphaU_imag = [];
alphaL_real = [];
alphaL_imag = [];
for k = 1:numIterations
    d = data(k);
    omega_real  = [omega_real, real(d.omega_F)];
    omega_imag  = [omega_imag, imag(d.omega_F)];
    L_real      = [L_real, real(d.L)];
    L_imag      = [L_imag, imag(d.L)];
    F_real      = [F_real, real(d.F)];
    F_imag      = [F_imag, imag(d.F)];
    alphaU_real = [alphaU_real, real(d.alpha_L_u)];
    alphaU_imag = [alphaU_imag, imag(d.alpha_L_u)];
    alphaL_real = [alphaL_real, real(d.alpha_L_l)];
    alphaL_imag = [alphaL_imag, imag(d.alpha_L_l)];
end
omega_xlim = [min([omega_real, L_real]), max([omega_real, L_real])];
omega_ylim = [min([omega_imag, L_imag]), max([omega_imag, L_imag])];
alpha_xlim = [min([F_real, alphaU_real, alphaL_real]), max([F_real, alphaU_real, alphaL_real])];
alpha_ylim = [min([F_imag, alphaU_imag, alphaL_imag]), max([F_imag, alphaU_imag, alphaL_imag])];
margin = 0.05;
omega_xlim = omega_xlim + margin * diff(omega_xlim) * [-1, 1];
omega_ylim = omega_ylim + margin * diff(omega_ylim) * [-1, 1];
alpha_xlim = alpha_xlim + margin * diff(alpha_xlim) * [-1, 1];
alpha_ylim = alpha_ylim + margin * diff(alpha_ylim) * [-1, 1];

make_contour_video(data, OUT_FILE, omega_xlim, omega_ylim, alpha_xlim, alpha_ylim, ...
                   ITER_START, ITER_STRIDE);

%% ======================================================================
%% Load JSON  --  THIS IS THE FIX
%% ======================================================================
function data = read_contour_integration_json(filename)

    % v6 onwards writes JSONL: one complete JSON object per line, appended.
    % Decode it line by line -- that also makes a run killed mid-write harmless,
    % because the partial final line is simply dropped.  The old single-array
    % .json files still load through the branch below, unchanged.
    [~, ~, ext] = fileparts(filename);
    if strcmpi(ext, '.jsonl')
        jsonData = read_jsonl(filename);
        if isempty(jsonData)
            error('read_contour:empty', '%s has no complete entries.', filename);
        end
    else

    raw = fileread(filename);

    try
        jsonData = jsondecode(raw);
    catch err
        error('read_contour:decode', ...
              ['jsondecode failed on %s (%s).\n' ...
               'Versions up to v5.1 rewrote this file in full on every save, ' ...
               'so a run killed mid-write leaves it truncated. Check that the ' ...
               'file ends with "}]" and fall back to a backup copy if not.'], ...
              filename, err.message);
    end

    % jsondecode gives a struct array only if every entry has the same field
    % names; iteration 1 has fewer, so v4.6/v4.7 come back as a cell array.
    % Normalise both to a cell array and index with {} below.
    if isstruct(jsonData)
        jsonData = num2cell(jsonData);
    elseif ~iscell(jsonData)
        error('read_contour:shape', ...
              'Top level of %s decoded as %s; expected an array of objects.', ...
              filename, class(jsonData));
    end

    end   % of the .json / .jsonl branch

    n = numel(jsonData);

    data = repmat(struct( ...
        'iteration',   [], ...
        'F',           [], ...
        'L',           [], ...
        'omega_F',     [], ...
        'alpha_L_u',   [], ...
        'alpha_L_l',   [], ...
        'd_branch',    [], ...
        'adapt_level', [], ...
        'n_L',         [], ...
        'zeta_alpha',  [], ...
        'exp_arg_max', []), n, 1);

    for k = 1:n
        e = jsonData{k};

        data(k).iteration   = get_scalar(e, 'iteration',   k);
        data(k).F           = reim_to_complex(get_field(e, 'F'));
        data(k).L           = reim_to_complex(get_field(e, 'L'));
        data(k).omega_F     = reim_to_complex(get_field(e, 'omega_F'));
        data(k).alpha_L_u   = reim_to_complex(get_field(e, 'alpha_L_u'));
        data(k).alpha_L_l   = reim_to_complex(get_field(e, 'alpha_L_l'));

        % scalars: absent on iteration 1, and null (-> []) whenever the Julia
        % side wrote a non-finite value through jsonnum().  NaN for both.
        data(k).d_branch    = get_scalar(e, 'd_branch',    NaN);
        data(k).adapt_level = get_scalar(e, 'adapt_level', NaN);
        data(k).n_L         = get_scalar(e, 'n_L',         NaN);
        data(k).zeta_alpha  = get_scalar(e, 'zeta_alpha',  NaN);
        data(k).exp_arg_max = get_scalar(e, 'exp_arg_max', NaN);
    end
end

function v = get_field(s, name)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = [];
    end
end

function x = get_scalar(s, name, fallback)
    x = fallback;
    if isstruct(s) && isfield(s, name)
        v = s.(name);
        if isnumeric(v) && isscalar(v)
            x = double(v);
        end
    end
end

%% ----------------------------------------------------------------------
%% ======================================================================
%% JSONL reader -- one object per line, preallocated, crash tolerant
%% ======================================================================
function c = read_jsonl(filename)

    fid = fopen(filename, 'r');
    if fid < 0
        error('read_contour:open', 'Could not open %s.', filename);
    end
    closer = onCleanup(@() fclose(fid));

    c = cell(4096, 1);
    n = 0;
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if isempty(strtrim(line))
            continue
        end
        try
            e = jsondecode(line);
        catch
            break      % partial final line from a killed run: stop cleanly
        end
        n = n + 1;
        if n > numel(c)
            c{2 * numel(c), 1} = [];
        end
        c{n} = e;
    end
    c = c(1:n);
end

%% ======================================================================
%% Convert the JSON {re, im} array to a MATLAB complex row vector
%% ----------------------------------------------------------------------
function cvec = reim_to_complex(reimArray)

    if isempty(reimArray)
        cvec = complex(zeros(1, 0));
        return
    end

    if isstruct(reimArray)
        % the normal case: an N-by-1 struct array with fields re and im
        cvec = reshape(complex([reimArray.re], [reimArray.im]), 1, []);
        return
    end

    if iscell(reimArray)
        % ragged fallback, in case a future writer emits a null component
        n = numel(reimArray);
        cvec = complex(zeros(1, n));
        for i = 1:n
            e = reimArray{i};
            re = 0; im = 0;
            if isstruct(e)
                if isfield(e, 're') && isscalar(e.re), re = e.re; end
                if isfield(e, 'im') && isscalar(e.im), im = e.im; end
            end
            cvec(i) = complex(re, im);
        end
        return
    end

    error('reim_to_complex:shape', ...
          'Expected a struct array of {re, im}, got %s.', class(reimArray));
end

%% ----------------------------------------------------------------------
%% Make contour video  (unchanged except: takes `data`, honours the stride,
%% and prints the diagnostics in the figure title)
%% ----------------------------------------------------------------------
function make_contour_video(data, outputfile, omega_xlim, omega_ylim, ...
                            alpha_xlim, alpha_ylim, iterStart, iterStride)

    numIterations = length(data);

    v = VideoWriter(outputfile, 'MPEG-4');
    v.FrameRate = 24;
    FrameRateLengthOnePlot = 2;
    open(v);
    cleanupObj = onCleanup(@() safe_close_video(v));

    figure('Position', [100, 100, 2000, 1000]);
    h1 = subplot(1, 2, 1);
    h2 = subplot(1, 2, 2);

    L      = data(iterStart).L;
    omega  = data(iterStart).omega_F;
    F      = data(iterStart).F;
    alphaU = data(iterStart).alpha_L_u;
    alphaL = data(iterStart).alpha_L_l;

    for k = iterStart:iterStride:numIterations
        d = data(k);
        % "\\zeta" so sprintf emits a literal \zeta for the TeX interpreter
        ttl = sprintf('iter %d    d = %.3e    \\zeta = %.2e    xarg = %.1f    level %g    N_L = %g', ...
                      d.iteration, d.d_branch, d.zeta_alpha, d.exp_arg_max, ...
                      d.adapt_level, d.n_L);

        % ---- Frame 1: update L ----
        L = d.L;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, ttl);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);

        % ---- Frame 2: update alpha_L_u and alpha_L_l ----
        alphaU = d.alpha_L_u;
        alphaL = d.alpha_L_l;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, ttl);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);

        % ---- Frame 3: update F ----
        F = d.F;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, ttl);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);

        % ---- Frame 4: update omega_F ----
        omega = d.omega_F;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, ttl);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);
    end

    close(v);
    disp('Video completed.');
end

function safe_close_video(v)
    % same as in plot_contour_deformation_v1.m
    try
        close(v);
    catch
        % it may already be closed; ignore
    end
end

function write_repeated_frames(v, figHandle, numRepeats)
    frame = getframe(figHandle);
    for i = 1:numRepeats
        writeVideo(v, frame);
    end
end

function plot_temporal(ax, L, omega_F, xlimits, ylimits, ttl)
    axes(ax); cla;
    hold on; grid on;
    plot(real(L), imag(L), 'b', 'DisplayName', 'L');
    plot(real(omega_F), imag(omega_F), 'r', 'DisplayName', '\omega_F');
    xlabel('\omega_r'); ylabel('\omega_i');
    title({'\omega-plane', ttl});
    legend show;
    axis normal;
    xlim(xlimits);
    ylim(ylimits);
end

function plot_spatial(ax, F, alphaU, alphaL, xlimits, ylimits)
    axes(ax); cla;
    hold on; grid on;
    plot(real(F), imag(F), 'b', 'DisplayName', 'F');
    plot(real(alphaU), imag(alphaU), 'r', 'DisplayName', '\alpha_L^u');
    plot(real(alphaL), imag(alphaL), 'g', 'DisplayName', '\alpha_L^l');
    xlabel('\alpha_r'); ylabel('\alpha_i');
    title('\alpha-plane');
    legend show;
    axis normal;
    xlim(xlimits);
    ylim(ylimits);
end
