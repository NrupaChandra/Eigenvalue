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
    for k = 1:numIterations-1
        entry = jsonData(k);
        data(k).iteration = entry.iteration;
        data(k).F         = reim_to_complex(entry.F);
        data(k).L         = reim_to_complex(entry.L);
        data(k).omega_F   = reim_to_complex(entry.omega_F);
        data(k).alpha_L_u = reim_to_complex(entry.alpha_L_u);
        data(k).alpha_L_l = reim_to_complex(entry.alpha_L_l);
    end
end

function cvec = reim_to_complex(reimArray)
    n = length(reimArray);
    cvec = complex(zeros(1, n));
    for i = 1:n
        cvec(i) = complex(reimArray(i).re, reimArray(i).im);
    end
end

data = read_contour_integration_json('contour_iteration.json');

%% Video

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

function make_contour_video(jsonfile, outputfile, omega_xlim, omega_ylim, alpha_xlim, alpha_ylim)
    data = read_contour_integration_json(jsonfile);
    numIterations = length(data);
    v = VideoWriter(outputfile, 'MPEG-4');
    v.FrameRate = 24;
    FrameRateLengthOnePlot = 2;
    open(v);
    figure('Position', [100, 100, 2000, 1000]);
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
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);  

        % ---- Frame 2: update alpha_L_u and alpha_L_l ----
        alphaU = d.alpha_L_u;
        alphaL = d.alpha_L_l;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);

        % ---- Frame 3: update F ----
        F = d.F;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);

        % ---- Frame 4: update omega_F ----
        omega = d.omega_F;
        plot_temporal(h1, L, omega, omega_xlim, omega_ylim);
        plot_spatial(h2, F, alphaU, alphaL, alpha_xlim, alpha_ylim);
        write_repeated_frames(v, gcf, FrameRateLengthOnePlot);
    end
    close(v);
    disp('Video completed.');
end

function write_repeated_frames(v, figHandle, numRepeats)
    frame = getframe(figHandle);
    for i = 1:numRepeats
        writeVideo(v, frame);
    end
end

function plot_temporal(ax, L, omega_F, xlimits, ylimits)
    axes(ax); cla;
    hold on; grid on;
    plot(real(L), imag(L), 'b', 'DisplayName', 'L');
    plot(real(omega_F), imag(omega_F), 'r', 'DisplayName', '\omega_F');
    xlabel('\omega_r'); ylabel('\omega_i');
    title('\omega-plane');
    legend show;
    axis equal;
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
    axis equal;
    xlim(xlimits);
    ylim(ylimits);
end

make_contour_video('contour_iteration.json', 'contour_video.mp4', omega_xlim, omega_ylim, alpha_xlim, alpha_ylim)