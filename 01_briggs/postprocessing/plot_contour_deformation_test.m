data = read_contour_integration_json('contour_iteration_algebraic_testv6.json');

%% Build fixed axis limits for the full video

numIterations = length(data);

all_omega = [];
all_alpha = [];

for k = 1:numIterations
    d = data(k);

    all_omega = [all_omega; d.omega_F(:); d.L(:)];
    all_alpha = [all_alpha; d.F(:); d.alpha_L_u(:); d.alpha_L_l(:)];
end

% Remove NaN / missing values
all_omega = all_omega(isfinite(real(all_omega)) & isfinite(imag(all_omega)));
all_alpha = all_alpha(isfinite(real(all_alpha)) & isfinite(imag(all_alpha)));

omega_xlim = [min(real(all_omega)), max(real(all_omega))];
omega_ylim = [min(imag(all_omega)), max(imag(all_omega))];

alpha_xlim = [min(real(all_alpha)), max(real(all_alpha))];
alpha_ylim = [min(imag(all_alpha)), max(imag(all_alpha))];

margin = 0.05;

omega_dx = omega_xlim(2) - omega_xlim(1);
omega_dy = omega_ylim(2) - omega_ylim(1);
alpha_dx = alpha_xlim(2) - alpha_xlim(1);
alpha_dy = alpha_ylim(2) - alpha_ylim(1);

% Avoid zero-width axes
if omega_dx == 0
    omega_dx = 1;
end
if omega_dy == 0
    omega_dy = 1;
end
if alpha_dx == 0
    alpha_dx = 1;
end
if alpha_dy == 0
    alpha_dy = 1;
end

omega_xlim = omega_xlim + margin * omega_dx * [-1, 1];
omega_ylim = omega_ylim + margin * omega_dy * [-1, 1];

alpha_xlim = alpha_xlim + margin * alpha_dx * [-1, 1];
alpha_ylim = alpha_ylim + margin * alpha_dy * [-1, 1];

%% Make video

make_contour_video( ...
    'contour_iteration_algebraic_testv6.json', ...
    'contour_video19_75frames.mp4', ...
    omega_xlim, omega_ylim, alpha_xlim, alpha_ylim ...
);

%% Load JSON

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

%% Convert JSON real-imag structure to complex vector

function cvec = reim_to_complex(reimArray)

    n = numel(reimArray);
    cvec = complex(nan(n,1), nan(n,1));

    for i = 1:n
        re = reimArray(i).re;
        im = reimArray(i).im;

        if isempty(re) || isempty(im)
            cvec(i) = NaN + 1i*NaN;
        else
            cvec(i) = complex(re, im);
        end
    end
end

%% Make contour video

function make_contour_video(jsonfile, outputfile, omega_xlim, omega_ylim, alpha_xlim, alpha_ylim)

    data = read_contour_integration_json(jsonfile);
    numIterations = length(data);

    v = VideoWriter(outputfile, 'MPEG-4');
    v.FrameRate = 75;

    FrameRateLengthOnePlot = 1;

    open(v);

    fig = figure('Position', [100, 100, 2000, 1000]);

    h1 = subplot(1, 2, 1);
    h2 = subplot(1, 2, 2);

    for k = 1:numIterations

        d = data(k);

        if k == 1
            L      = d.L;
            omega  = d.omega_F;
            F      = d.F;
            alphaU = d.alpha_L_u;
            alphaL = d.alpha_L_l;
        end

        % ---- Frame 1: update L ----
        L = d.L;

        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, k);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim, k);

        write_repeated_frames(v, fig, FrameRateLengthOnePlot);

        % ---- Frame 2: update alpha_L_u and alpha_L_l ----
        alphaU = d.alpha_L_u;
        alphaL = d.alpha_L_l;

        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, k);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim, k);

        write_repeated_frames(v, fig, FrameRateLengthOnePlot);

        % ---- Frame 3: update F ----
        F = d.F;

        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, k);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim, k);

        write_repeated_frames(v, fig, FrameRateLengthOnePlot);

        % ---- Frame 4: update omega_F ----
        omega = d.omega_F;

        plot_temporal(h1, L, omega, omega_xlim, omega_ylim, k);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim, k);

        write_repeated_frames(v, fig, FrameRateLengthOnePlot);
    end

    close(v);
    close(fig);

    disp(['Video completed: ', outputfile]);
end

%% Write repeated frames

function write_repeated_frames(v, figHandle, numRepeats)

    drawnow;
    frame = getframe(figHandle);

    for i = 1:numRepeats
        writeVideo(v, frame);
    end
end

%% Plot omega plane

function plot_temporal(ax, L, omega_F, xlimits, ylimits, k)

    axes(ax);
    cla;
    hold on;
    grid on;

    valid_L = isfinite(real(L)) & isfinite(imag(L));
    valid_omega = isfinite(real(omega_F)) & isfinite(imag(omega_F));

    plot(real(L(valid_L)), imag(L(valid_L)), ...
        'b.-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 8, ...
        'DisplayName', 'L');

    plot(real(omega_F(valid_omega)), imag(omega_F(valid_omega)), ...
        'r.-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 8, ...
        'DisplayName', '\omega_F');

    % Exact omega pinch marker for this algebraic test
    omega_pinch_exact = -8.800000000000e-01 +8.000000000000e-01i;

    plot(real(omega_pinch_exact), imag(omega_pinch_exact), ...
        'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 7, ...
        'DisplayName', 'exact \omega pinch');

    xlabel('\omega_r');
    ylabel('\omega_i');
    title(sprintf('\\omega-plane, iteration %d', k));

    legend('Location', 'best');

    axis normal;
    xlim(xlimits);
    ylim(ylimits);
end

%% Plot alpha plane

function plot_spatial(ax, F, alphaU, alphaL, xlimits, ylimits, k)

    axes(ax);
    cla;
    hold on;
    grid on;

    valid_F = isfinite(real(F)) & isfinite(imag(F));
    valid_U = isfinite(real(alphaU)) & isfinite(imag(alphaU));
    valid_L = isfinite(real(alphaL)) & isfinite(imag(alphaL));

    plot(real(F(valid_F)), imag(F(valid_F)), ...
        'b.-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 8, ...
        'DisplayName', 'F');

    plot(real(alphaU(valid_U)), imag(alphaU(valid_U)), ...
        'r.-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 8, ...
        'DisplayName', '\alpha_L^u');

    plot(real(alphaL(valid_L)), imag(alphaL(valid_L)), ...
        'g.-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 8, ...
        'DisplayName', '\alpha_L^l');

    % Exact alpha pinch marker for this algebraic test
    alpha_pinch_exact =-2.000000000000e-01 -4.000000000000e-01i;

    plot(real(alpha_pinch_exact), imag(alpha_pinch_exact), ...
        'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 7, ...
        'DisplayName', 'exact \alpha pinch');

    xlabel('\alpha_r');
    ylabel('\alpha_i');
    %title(sprintf('\alpha-plane, iteration %d', k));
    hTitle = title(ax, sprintf('\\alpha-plane, iteration %d', k));
    hTitle.Units = 'normalized';
    hTitle.Position = [0.5, 1.02, 0];
    hTitle.HorizontalAlignment = 'center';

   % legend('Location', 'best');
    legend(ax, 'Location', 'northeast');
    axis normal;
    xlim(xlimits);
    ylim(ylimits);
end