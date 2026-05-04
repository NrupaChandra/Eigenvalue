data = read_contour_integration_json('contour_iteration.json');

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

    if isempty(d.L) || isempty(d.F) || isempty(d.omega_F) || ...
       isempty(d.alpha_L_u) || isempty(d.alpha_L_l)
        fprintf('Skipping empty entry while computing limits: k = %d\n', k);
        continue
    end

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

% make_contour_video( ...   
%     'contour_iteration.json', ...
%     'contour_video_partial_saved.mp4', ...
%     omega_xlim, omega_ylim, ...
%     alpha_xlim, alpha_ylim ...
% );
make_contour_video( ...
    'contour_iteration.json', ...
    'contour_video_partial_saved.avi', ...
    omega_xlim, omega_ylim, ...
    alpha_xlim, alpha_ylim ...
);
%% ------------------------------------------------------------------------
%% Load JSON
%% ------------------------------------------------------------------------

function data = read_contour_integration_json(filename)

    raw = fileread(filename);
    jsonData = jsondecode(raw);

    numIterations = length(jsonData);

    data(numIterations) = struct( ...
        'iteration', [], ...
        'F', [], ...
        'L', [], ...
        'omega_F', [], ...
        'alpha_L_u', [], ...
        'alpha_L_l', [] ...
    );

    for k = 1:numIterations
        entry = jsonData(k);

        data(k).iteration = entry.iteration;
        data(k).F         = reim_to_complex(entry.F);
        data(k).L         = reim_to_complex(entry.L);
        data(k).omega_F   = reim_to_complex(entry.omega_F);
        data(k).alpha_L_u = reim_to_complex(entry.alpha_L_u);
        data(k).alpha_L_l = reim_to_complex(entry.alpha_L_l);
    end

end

%% ------------------------------------------------------------------------
%% Convert JSON real-imag format to MATLAB complex vector
%% ------------------------------------------------------------------------

function cvec = reim_to_complex(reimArray)

    n = length(reimArray);
    cvec = complex(zeros(1, n));

    for i = 1:n
        cvec(i) = complex(reimArray(i).re, reimArray(i).im);
    end

end

%% ------------------------------------------------------------------------
%% Make contour video
%% ------------------------------------------------------------------------

function make_contour_video(jsonfile, outputfile, omega_xlim, omega_ylim, alpha_xlim, alpha_ylim)

    data = read_contour_integration_json(jsonfile);
    numIterations = length(data);

    % v = VideoWriter(outputfile, 'MPEG-4');
    % v.FrameRate = 24;
    % v.Quality = 80;
    v = VideoWriter(outputfile, 'Motion JPEG AVI');
    v.FrameRate = 24;
    v.Quality = 75;

    FrameRateLengthOnePlot = 1;

    open(v);

    % This tries to finalize the mp4 if an error happens.
    cleanupObj = onCleanup(@() safe_close_video(v));

    fig = figure('Position', [100, 100, 1200, 600]);

    h1 = subplot(1, 2, 1); 
    h2 = subplot(1, 2, 2); 

    L = [];
    omega = [];
    F = [];
    alphaU = [];
    alphaL = [];

    for k = 1:numIterations

        try
            d = data(k);

            if isempty(d.L) || isempty(d.F) || isempty(d.omega_F) || ...
               isempty(d.alpha_L_u) || isempty(d.alpha_L_l)
                fprintf('Skipping empty iteration k = %d\n', k);
                continue
            end

            if isempty(L)
                L      = d.L;
                omega  = d.omega_F;
                F      = d.F;
                alphaU = d.alpha_L_u;
                alphaL = d.alpha_L_l;
            end

            % ---- Frame 1: update L ----
            L = d.L;
            plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
            plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
            write_repeated_frames(v, fig, FrameRateLengthOnePlot);  

            % ---- Frame 2: update alpha_L_u and alpha_L_l ----
            alphaU = d.alpha_L_u;
            alphaL = d.alpha_L_l;
            plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
            plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
            write_repeated_frames(v, fig, FrameRateLengthOnePlot);

            % ---- Frame 3: update F ----
            F = d.F;
            plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
            plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
            write_repeated_frames(v, fig, FrameRateLengthOnePlot);

            % ---- Frame 4: update omega_F ----
            omega = d.omega_F;
            plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
            plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
            write_repeated_frames(v, fig, FrameRateLengthOnePlot);

            fprintf('Saved iteration %d / %d\n', k, numIterations);

        catch ME
            fprintf('\nVideo creation crashed at iteration %d.\n', k);
            fprintf('Reason: %s\n', ME.message);
            fprintf('Trying to close/finalize video now...\n');

            safe_close_video(v);

            if isvalid(fig)
                close(fig);
            end

            return
        end
    end

    safe_close_video(v);

    if isvalid(fig)
        close(fig);
    end

    disp('Video completed.');

end

%% ------------------------------------------------------------------------
%% Safely close video
%% ------------------------------------------------------------------------

function safe_close_video(v)

    try
        close(v);
        disp('Video safely closed/finalized.');
    catch
        % It may already be closed. Ignore.
    end

end

%% ------------------------------------------------------------------------
%% Write repeated frames
%% ------------------------------------------------------------------------

function write_repeated_frames(v, figHandle, numRepeats)

    drawnow limitrate;
    frame = getframe(figHandle);

    for i = 1:numRepeats
        writeVideo(v, frame);
    end

end

%% ------------------------------------------------------------------------
%% Plot temporal omega-plane
%% ------------------------------------------------------------------------

function plot_temporal(ax, L, omega_F, xlimits, ylimits)

    axes(ax); 
    cla(ax);

    hold(ax, 'on'); 
    grid(ax, 'on');

    plot(ax, real(L), imag(L), 'b', 'DisplayName', 'L');
    plot(ax, real(omega_F), imag(omega_F), 'r', 'DisplayName', '\omega_F');

    xlabel(ax, '\omega_r'); 
    ylabel(ax, '\omega_i');

    title(ax, '\omega-plane');

    legend(ax, 'show');

    axis(ax, 'normal');

    xlim(ax, xlimits);
    ylim(ax, ylimits);

end

%% ------------------------------------------------------------------------
%% Plot spatial alpha-plane
%% ------------------------------------------------------------------------

function plot_spatial(ax, F, alphaU, alphaL, xlimits, ylimits)

    axes(ax); 
    cla(ax);

    hold(ax, 'on'); 
    grid(ax, 'on');

    plot(ax, real(F), imag(F), 'b', 'DisplayName', 'F');
    plot(ax, real(alphaU), imag(alphaU), 'r', 'DisplayName', '\alpha_L^u');
    plot(ax, real(alphaL), imag(alphaL), 'g', 'DisplayName', '\alpha_L^l');

    xlabel(ax, '\alpha_r'); 
    ylabel(ax, '\alpha_i');

    title(ax, '\alpha-plane');

    legend(ax, 'show');

    axis(ax, 'normal');

    xlim(ax, xlimits);
    ylim(ax, ylimits);

end