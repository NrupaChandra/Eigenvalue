clear; clc;
close all;

% Read data once
filename = 'combined_data_standardized.json'; % JSON file containing your data
jsonData = jsondecode(fileread(filename));  % Read and decode JSON data

keepGoing = true;  % Control flag for the menu loop

while keepGoing
    % Display available analysis options including a Quit option
    fprintf('\nAvailable Analysis Options:\n');
    fprintf('  1: Plot eigenvalues for a specific Re and omega\n');
    fprintf('  2: Compare ratios (omega/eigenvalue) with group velocity (with optional filtering by Re)\n');
    fprintf('  3: Plot unstable eigenvalues for a specified range of Re\n');
    fprintf('  4: Plot unstable eigenvalues for a specified omega range (for a given Re)\n');
    fprintf('  5: Combined plot (target Re & omega with unstable modes over a Re range)\n');
    fprintf('  6: Plot unstable eigenvalues for all Re for a given range of omega (real part)\n');
    fprintf('  7: Combined plot (target Re & omega with unstable modes over a range of omega)\n');
    fprintf('  8: Quit\n');
    
    % Get user input
    option = input('\nEnter your choice (1-8): ');
    
    % Validate the input (if not 1 to 8, notify the user and display menu again)
    if isempty(option) || ~ismember(option, 1:8)
        fprintf('Invalid option. Please try again.\n');
        continue;  % Go back to the beginning of the loop
    end
    
    % Quit option
    if option == 8
        fprintf('Quitting the program. Goodbye!\n');
        break;  % Exit the while loop
    end

    % Main Switch to the corresponding function based on user selection
    switch option
        case 1
            myRe    = input('Enter target Reynolds number (e.g., 1600): ');
            myOmega = input('Enter target omega (as complex, e.g., 0.02-0.3i): ');
            plotAlphaReIm(jsonData, myRe, myOmega);
            
        case 2
            applyFilter = input('Do you want to filter by a specific Reynolds number? (1 for Yes, 0 for No): ');
            if applyFilter
                targetRe = input('Enter target Reynolds number for filtering: ');
                compareRatiosWithGroupVelocityCombined(jsonData, true, targetRe);
            else
                compareRatiosWithGroupVelocityCombined(jsonData, false, 0);
            end
            
        case 3
            minRe = input('Enter minimum Reynolds number: ');
            maxRe = input('Enter maximum Reynolds number: ');
            plotUnstableEigenvalues(jsonData, minRe, maxRe);
            
        case 4
            myRe       = input('Enter target Reynolds number: ');
            minOmegaRe = input('Enter minimum omega (real part): ');
            maxOmegaRe = input('Enter maximum omega (real part): ');
            plotUnstableEigenvaluesForOmega(jsonData, myRe, minOmegaRe, maxOmegaRe);
            
        case 5
            % Combined plot: target case eigenvalues with unstable modes over a Re range
            myRe    = input('Enter target Reynolds number (e.g., 1600): ');
            myOmega = input('Enter target omega (as complex, e.g., 0.02+0.1i): ');
            minRe   = input('Enter minimum Reynolds number for unstable modes: ');
            maxRe   = input('Enter maximum Reynolds number for unstable modes: ');
            plotCombinedEigenvalues(jsonData, myRe, myOmega, minRe, maxRe);
            
        case 6
            % Plot unstable eigenvalues for all Re for a given omega range (real part)
            minOmegaRe = input('Enter minimum omega (real part): ');
            maxOmegaRe = input('Enter maximum omega (real part): ');
            plotUnstableEigenvaluesAllReForOmegaRange(jsonData, minOmegaRe, maxOmegaRe);
            
        case 7
            % Combined Option 7: Merge Option 1 and Option 6 into one plot.
            myRe    = input('Enter target Reynolds number (e.g., 1600): ');
            myOmega = input('Enter target omega (as complex, e.g., 0.02+0.1i): ');
            minOmegaRe = input('Enter minimum omega (real part) for unstable eigenvalue plot: ');
            maxOmegaRe = input('Enter maximum omega (real part) for unstable eigenvalue plot: ');
            plotCombinedOption7(jsonData, myRe, myOmega, minOmegaRe, maxOmegaRe);
            
        otherwise
            fprintf('Invalid option. Exiting script.\n');
            break;
    end
    
    % Ask if the user wants to perform another analysis
    cont_choice = input('\nWould you like to choose another option? (y/n): ', 's');
    if lower(cont_choice) ~= 'y'
        fprintf('Quitting the program. Goodbye!\n');
        break;
    end
end

%% Function Definitions
function plotAlphaReIm(jsonData, myRe, myOmega)
    % Compute the implied alpha based on the target Reynolds number.
    implied_phi = 1432.0 / myRe;
    implied_alpha_2d_r = 2 * pi / 50 * implied_phi;
    implied_alpha_2d_i = 2 * pi / 25 * 1 / (sqrt(1/(implied_phi^2) - 1));
    tol = 1e-6;
    
    % Loop through JSON data to find matching Re and omega.
    for i = 1:numel(jsonData)
        data = jsonData(i);
        % Extract omega (handles both numeric and structured formats)
        if isnumeric(data.omega)
            omega_real = data.omega;
            omega_imag = 0;
        elseif isstruct(data.omega)
            omega_real = data.omega.re;
            omega_imag = data.omega.im;
        else
            error('Unknown format for omega.');
        end
        dataOmega = omega_real + 1i * omega_imag;
        
        if (data.reynolds == myRe) && (abs(real(dataOmega)-real(myOmega)) < tol) && (abs(imag(dataOmega)-imag(myOmega)) < tol)
            ev_iterated1 = data.ev_iterated;
            ev_iterated_re = [];
            ev_iterated_im = [];
            value_instab_2d = [];
            value_instab_3d = [];
            for j = 1:numel(ev_iterated1)
                ev = ev_iterated1(j);
                if isfield(ev.value, 're') && ~isempty(ev.value.re)
                    ev_iterated_re(j) = ev.value.re;
                else
                    ev_iterated_re(j) = NaN;
                end
                if isfield(ev.value, 'im') && ~isempty(ev.value.im)
                    ev_iterated_im(j) = ev.value.im;
                else
                    ev_iterated_im(j) = NaN;
                end
                if isfield(ev, 'value_instab_2d') && ~isempty(ev.value_instab_2d)
                    value_instab_2d(j) = ev.value_instab_2d(1);
                else
                    value_instab_2d(j) = NaN;
                end
                if isfield(ev, 'value_instab_3d') && ~isempty(ev.value_instab_3d)
                    value_instab_3d(j) = ev.value_instab_3d(1);
                else
                    value_instab_3d(j) = NaN;
                end
            end
            
            % Create scatter plot for eigenvalues.
            figure;
            hold on; grid on;
            plot(ev_iterated_re, ev_iterated_im, 'kd', 'MarkerSize', 4, 'DisplayName', 'Eigenvalues');
            for j = 1:numel(ev_iterated1)
                if value_instab_2d(j) < 0
                    plot(ev_iterated_re(j), ev_iterated_im(j), 'bd', 'MarkerSize', 8, 'DisplayName', 'Unstable 2D');
                end
                if value_instab_3d(j) < 0
                    plot(ev_iterated_re(j), ev_iterated_im(j), 'md', 'MarkerSize', 12, 'DisplayName', 'Unstable 3D');
                end
            end
            plot(implied_alpha_2d_r, implied_alpha_2d_i, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Implied alpha');
            xlabel('\alpha_{2D, real}');
            ylabel('\alpha_{2D, imag}');
            title(['Eigenvalues for Re = ', num2str(myRe), ', \omega = ', num2str(myOmega)]);
            hold off;
        end
    end
end

function compareRatiosWithGroupVelocityCombined(jsonData, applyFilter, targetRe)
    % Compute the ratio (omega/eigenvalue) and plot its relationship with group velocity.
    group_velocity_reals = [];
    group_velocity_imags = [];
    ratio_values_real = [];
    ratio_values_imaginary = [];
    
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if applyFilter && (data.reynolds ~= targetRe)
            continue;
        end
        
        if isnumeric(data.omega)
            omega_real = data.omega;
            omega_imag = 0;
        elseif isstruct(data.omega)
            omega_real = data.omega.re;
            omega_imag = data.omega.im;
        else
            error('Unknown format for omega.');
        end
        
        for j = 1:numel(data.ev_iterated)
            ev = data.ev_iterated(j).value;
            if isfield(ev, 're') && ~isempty(ev.re)
                ev_real = ev.re;
                ev_imag = ev.im;
                group_velocity_real = data.ev_iterated(j).group_velocity.re;
                group_velocity_imag = data.ev_iterated(j).group_velocity.im;
                ratio_omega_Re = omega_real / ev_real;
                if ~isempty(omega_imag) && ~isempty(ev_imag)
                    ratio_omega_Im = omega_imag / ev_imag;
                else
                    ratio_omega_Im = NaN;
                end
                group_velocity_reals(end+1,1) = group_velocity_real;
                group_velocity_imags(end+1,1) = group_velocity_imag;
                ratio_values_real(end+1,1) = ratio_omega_Re;
                ratio_values_imaginary(end+1,1) = ratio_omega_Im;
            end
        end
    end
    
    figure;
    subplot(1,2,1);
    scatter(ratio_values_real, group_velocity_reals, 'filled');
    xlabel('omega/Eigenvalue (Real)');
    ylabel('Group Velocity (Real)');
    title('Real Parts Comparison');
    grid on;
    subplot(1,2,2);
    scatter(ratio_values_imaginary, group_velocity_imags, 'filled');
    xlabel('omega/Eigenvalue (Imag)');
    ylabel('Group Velocity (Imag)');
    title('Imaginary Parts Comparison');
    grid on;
end

function plotUnstableEigenvalues(jsonData, minRe, maxRe)
    % Plot unstable eigenvalues for data with Reynolds numbers between minRe and maxRe.
    unstable_2d_alpha_r = [];
    unstable_2d_alpha_i = [];
    unstable_3d_alpha_r = [];
    unstable_3d_alpha_i = [];
    
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if data.reynolds < minRe || data.reynolds > maxRe
            continue;
        end
        for j = 1:numel(data.ev_iterated)
            ev = data.ev_iterated(j);
            if isfield(ev.value, 're') && isfield(ev.value, 'im')
                alpha_r = ev.value.re;
                alpha_i = ev.value.im;
            else
                alpha_r = NaN; alpha_i = NaN;
            end
            if isfield(ev, 'value_instab_2d') && any(ev.value_instab_2d < 0)
                unstable_2d_alpha_r(end+1) = alpha_r;
                unstable_2d_alpha_i(end+1) = alpha_i;
                if isnumeric(data.omega)
                    omega_real = data.omega; omega_imag = 0;
                elseif isstruct(data.omega)
                    omega_real = data.omega.re; omega_imag = data.omega.im;
                else
                    error('Unknown format for omega.');
                end
                fprintf('2D Unstable Eigenvalue at Re: %g, omega: %g+%gi\n', data.reynolds, omega_real, omega_imag);
            end
            if isfield(ev, 'value_instab_3d') && any(ev.value_instab_3d < 0)
                unstable_3d_alpha_r(end+1) = alpha_r;
                unstable_3d_alpha_i(end+1) = alpha_i;
                if isnumeric(data.omega)
                    omega_real = data.omega; omega_imag = 0;
                elseif isstruct(data.omega)
                    omega_real = data.omega.re; omega_imag = data.omega.im;
                else
                    error('Unknown format for omega.');
                end
                fprintf('3D Unstable Eigenvalue at Re: %g, omega: %g+%gi\n', data.reynolds, omega_real, omega_imag);
            end
        end
    end
    
    figure;
    hold on; grid on;
    if ~isempty(unstable_2d_alpha_r)
        plot(unstable_2d_alpha_r, unstable_2d_alpha_i, 'bd', 'MarkerSize', 6, 'DisplayName', 'Unstable 2D');
    end
    if ~isempty(unstable_3d_alpha_r)
        plot(unstable_3d_alpha_r, unstable_3d_alpha_i, 'md', 'MarkerSize', 8, 'DisplayName', 'Unstable 3D');
    end
    xlabel('\alpha_{real}');
    ylabel('\alpha_{imag}');
    title(sprintf('Unstable Eigenvalues for Re range [%g, %g]', minRe, maxRe));
    hold off;
end

function plotUnstableEigenvaluesForOmega(jsonData, myRe, minOmegaRe, maxOmegaRe)
    % Plot unstable eigenvalues for a given Reynolds number (myRe) where the
    % real part of omega is within [minOmegaRe, maxOmegaRe].
    unstable_2d_alpha_r = [];
    unstable_2d_alpha_i = [];
    unstable_3d_alpha_r = [];
    unstable_3d_alpha_i = [];
    omega_real_values = [];
    omega_imag_values = [];
    
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if data.reynolds ~= myRe
            continue;
        end
        if isnumeric(data.omega)
            omega_real = data.omega;
            omega_imag = 0;
        elseif isstruct(data.omega)
            omega_real = data.omega.re;
            omega_imag = data.omega.im;
        else
            error('Unknown format for omega.');
        end
        if omega_real < minOmegaRe || omega_real > maxOmegaRe
            continue;
        end
        for j = 1:numel(data.ev_iterated)
            ev = data.ev_iterated(j);
            if isfield(ev.value, 're') && isfield(ev.value, 'im')
                alpha_r = ev.value.re;
                alpha_i = ev.value.im;
            else
                alpha_r = NaN; alpha_i = NaN;
            end
            if isfield(ev, 'value_instab_2d') && any(ev.value_instab_2d < 0)
                unstable_2d_alpha_r(end+1) = alpha_r;
                unstable_2d_alpha_i(end+1) = alpha_i;
                omega_real_values(end+1) = omega_real;
                omega_imag_values(end+1) = omega_imag;
                fprintf('2D Unstable at Re %g, omega: %g+%gi\n', myRe, omega_real, omega_imag);
            end
            if isfield(ev, 'value_instab_3d') && any(ev.value_instab_3d < 0)
                unstable_3d_alpha_r(end+1) = alpha_r;
                unstable_3d_alpha_i(end+1) = alpha_i;
                omega_real_values(end+1) = omega_real;
                omega_imag_values(end+1) = omega_imag;
                fprintf('3D Unstable at Re %g, omega: %g+%gi\n', myRe, omega_real, omega_imag);
            end
        end
    end
    
    figure;
    hold on; grid on;
    for k = 1:length(unstable_2d_alpha_r)
        plot(unstable_2d_alpha_r(k), unstable_2d_alpha_i(k), 'bd', 'MarkerSize', 6, ...
             'DisplayName', sprintf('omega = %g+%gi, 2D', omega_real_values(k), omega_imag_values(k)));
    end
    for k = 1:length(unstable_3d_alpha_r)
        plot(unstable_3d_alpha_r(k), unstable_3d_alpha_i(k), 'md', 'MarkerSize', 8, ...
             'DisplayName', sprintf('omega = %g+%gi, 3D', omega_real_values(k), omega_imag_values(k)));
    end
    xlabel('\alpha_{real}');
    ylabel('\alpha_{imag}');
    title(sprintf('Unstable Eigenvalues for Re = %g and omega_{real} in [%g, %g]', myRe, minOmegaRe, maxOmegaRe));
    hold off;
end

function plotCombinedEigenvalues(jsonData, myRe, myOmega, minRe, maxRe)
    % Define constants for implied alpha calculation
    implied_phi = 1400.0 / myRe;
    implied_alpha_2d_r = 2 * pi / 50 * implied_phi;
    implied_alpha_2d_i = 2 * pi / 25 * 1 / (sqrt(1/implied_phi^2 - 1));
    tol = 1e-6;

    % Initialize arrays to store unstable eigenvalues
    ev_iterated_re = [];
    ev_iterated_im = [];
    unstable_2d_alpha_r = [];
    unstable_2d_alpha_i = [];
    unstable_3d_alpha_r = [];
    unstable_3d_alpha_i = [];
    reynolds_values = [];
    
    % Get unique Reynolds numbers for the color map
    unique_reynolds = unique(arrayfun(@(x) x.reynolds, jsonData));
    cmap = colormap(jet(numel(unique_reynolds)));

    % Loop through each data entry in JSON
    for i = 1:numel(jsonData)
        data = jsonData(i);
        
        % Recreate complex omega from JSON data
        if isnumeric(data.omega)
            omega_real = data.omega;
            omega_imag = 0;
        elseif isstruct(data.omega)
            omega_real = data.omega.re;
            omega_imag = data.omega.im;
        else
            error('Unknown format for omega.');
        end
        dataOmega = omega_real + 1i * omega_imag;
        
        % Condition for specific Re and omega
        if (data.reynolds == myRe) && (abs(real(dataOmega) - real(myOmega)) < tol) && (abs(imag(dataOmega) - imag(myOmega)) < tol)
            for j = 1:numel(data.ev_iterated)
                ev = data.ev_iterated(j);
                if isfield(ev.value, 're') && isfield(ev.value, 'im')
                    if isscalar(ev.value.re) && isscalar(ev.value.im)
                        ev_iterated_re(end+1) = ev.value.re;
                        ev_iterated_im(end+1) = ev.value.im;
                    else
                        warning('Skipping non-scalar eigenvalue at index %d', j);
                    end
                end
            end
        end

        % Condition for unstable eigenvalues within specified Re range
        if data.reynolds >= minRe && data.reynolds <= maxRe
            for j = 1:numel(data.ev_iterated)
                ev = data.ev_iterated(j);
                if isfield(ev.value, 're') && isfield(ev.value, 'im')
                    if isfield(ev, 'value_instab_2d') && isscalar(ev.value_instab_2d) && (ev.value_instab_2d < 0)
                        unstable_2d_alpha_r(end+1) = ev.value.re;
                        unstable_2d_alpha_i(end+1) = ev.value.im;
                        reynolds_values(end+1) = data.reynolds;
                    end
                    if isfield(ev, 'value_instab_3d') && isscalar(ev.value_instab_3d) && (ev.value_instab_3d < 0)
                        unstable_3d_alpha_r(end+1) = ev.value.re;
                        unstable_3d_alpha_i(end+1) = ev.value.im;
                        reynolds_values(end+1) = data.reynolds;
                    end
                end
            end
        end
    end

    % Create combined plot
    figure;
    hold on;
    grid on;

    plot(ev_iterated_re, ev_iterated_im, 'd', 'MarkerSize', 4, 'Color', 'k', 'DisplayName', 'Reiterated eigenvalues');
    plot(implied_alpha_2d_r, implied_alpha_2d_i, 'x', 'MarkerSize', 4, 'Color', 'k', 'DisplayName', 'Implied alpha');

    % Plot unstable 2D and 3D eigenvalues, color-coded by Re
    for k = 1:numel(unstable_2d_alpha_r)
        color_idx = find(unique_reynolds == reynolds_values(k), 1);
        plot(unstable_2d_alpha_r(k), unstable_2d_alpha_i(k), 'd', 'MarkerSize', 6, 'Color', cmap(color_idx,:), 'DisplayName', ['Re = ' num2str(reynolds_values(k)) ', 2D']);
    end

    for k = 1:numel(unstable_3d_alpha_r)
        color_idx = find(unique_reynolds == reynolds_values(k), 1);
        plot(unstable_3d_alpha_r(k), unstable_3d_alpha_i(k), 'd', 'MarkerSize', 8, 'Color', cmap(color_idx,:), 'DisplayName', ['Re = ' num2str(reynolds_values(k)) ', 3D']);
    end

    xlabel('\alpha_{r}');
    ylabel('\alpha_{i}');
    title(['Combined Eigenvalues for Re = ', num2str(myRe), ', \omega = ', num2str(real(myOmega)), ' + ', num2str(imag(myOmega)), 'i']);
    colormap(cmap);
    colorbar('Ticks', linspace(0, 1, numel(unique_reynolds)), 'TickLabels', unique_reynolds);
    hold off;
end

function plotUnstableEigenvaluesAllReForOmegaRange(jsonData, minOmegaRe, maxOmegaRe)
    % Plot unstable eigenvalues for all Re where the real part of omega is in [minOmegaRe, maxOmegaRe].
    unstable_2d_alpha_r = [];
    unstable_2d_alpha_i = [];
    unstable_3d_alpha_r = [];
    unstable_3d_alpha_i = [];
    reynolds_2d = [];
    reynolds_3d = [];
    
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if isnumeric(data.omega)
            omega_real = data.omega;
            omega_imag = 0;
        elseif isstruct(data.omega)
            omega_real = data.omega.re;
            omega_imag = data.omega.im;
        else
            error('Unknown format for omega.');
        end
        
        if omega_real < minOmegaRe || omega_real > maxOmegaRe
            continue;
        end
        
        for j = 1:numel(data.ev_iterated)
            ev = data.ev_iterated(j);
            if isfield(ev.value, 're') && isfield(ev.value, 'im')
                alpha_r = ev.value.re;
                alpha_i = ev.value.im;
            else
                alpha_r = NaN; 
                alpha_i = NaN;
            end
            
            if isfield(ev, 'value_instab_2d') && any(ev.value_instab_2d < 0)
                unstable_2d_alpha_r(end+1) = alpha_r;
                unstable_2d_alpha_i(end+1) = alpha_i;
                reynolds_2d(end+1) = data.reynolds;
                fprintf('2D Unstable at Re %g, omega: %g+%gi\n', data.reynolds, omega_real, omega_imag);
            end
            
            if isfield(ev, 'value_instab_3d') && any(ev.value_instab_3d < 0)
                unstable_3d_alpha_r(end+1) = alpha_r;
                unstable_3d_alpha_i(end+1) = alpha_i;
                reynolds_3d(end+1) = data.reynolds;
                fprintf('3D Unstable at Re %g, omega: %g+%gi\n', data.reynolds, omega_real, omega_imag);
            end
        end
    end

    all_re = unique([reynolds_2d, reynolds_3d]);
    cmap = colormap(jet(numel(all_re)));
    
    figure;
    hold on; grid on;
    
    for k = 1:length(unstable_2d_alpha_r)
        color_idx = find(all_re == reynolds_2d(k), 1);
        plot(unstable_2d_alpha_r(k), unstable_2d_alpha_i(k), 'bd', 'MarkerSize', 6, 'Color', cmap(color_idx,:),...
             'DisplayName', sprintf('Re = %g, 2D', reynolds_2d(k)));
    end
    
    for k = 1:length(unstable_3d_alpha_r)
        color_idx = find(all_re == reynolds_3d(k), 1);
        plot(unstable_3d_alpha_r(k), unstable_3d_alpha_i(k), 'md', 'MarkerSize', 8, 'Color', cmap(color_idx,:),...
             'DisplayName', sprintf('Re = %g, 3D', reynolds_3d(k)));
    end
    
    xlabel('\alpha_{real}');
    ylabel('\alpha_{imag}');
    title(sprintf('Unstable Eigenvalues for All Re with omega_{real} in [%g, %g]', minOmegaRe, maxOmegaRe));
    colormap(cmap);
    colorbar('Ticks', linspace(0, 1, numel(all_re)), 'TickLabels', all_re);
    hold off;
end

function plotCombinedOption7(jsonData, myRe, myOmega, minOmegaRe, maxOmegaRe)
    % Combined Option 7: Plot specific eigenvalues (as in Option 1) along with
    % unstable eigenvalues for all Re with omega real part in the given range (as in Option 6)
    tol = 1e-6;
    
    eigen_re = [];
    eigen_im = [];
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if isnumeric(data.omega)
            curOmega = data.omega + 0i;
        elseif isstruct(data.omega)
            curOmega = data.omega.re + 1i*data.omega.im;
        else
            error('Unknown format for omega.');
        end
        
        if (data.reynolds == myRe) && (abs(real(curOmega)-real(myOmega)) < tol) && (abs(imag(curOmega)-imag(myOmega)) < tol)
            for j = 1:numel(data.ev_iterated)
                ev = data.ev_iterated(j);
                if isfield(ev.value, 're') && isfield(ev.value, 'im')
                    eigen_re(end+1) = ev.value.re;
                    eigen_im(end+1) = ev.value.im;
                end
            end
        end
    end
    
    implied_phi = 1432.0 / myRe;
    implied_alpha_2d_r = 2 * pi / 50 * implied_phi;
    implied_alpha_2d_i = 2 * pi / 25 * 1 / (sqrt(1/(implied_phi^2) - 1));
    
    unstable2d_re = [];
    unstable2d_im = [];
    unstable3d_re = [];
    unstable3d_im = [];
    reynolds_2d = [];
    reynolds_3d = [];
    
    for i = 1:numel(jsonData)
        data = jsonData(i);
        if isnumeric(data.omega)
            curOmega = data.omega + 0i;
        elseif isstruct(data.omega)
            curOmega = data.omega.re + 1i*data.omega.im;
        else
            error('Unknown format for omega.');
        end
        
        if real(curOmega) < minOmegaRe || real(curOmega) > maxOmegaRe
            continue;
        end
        
        for j = 1:numel(data.ev_iterated)
            ev = data.ev_iterated(j);
            if isfield(ev.value, 're') && isfield(ev.value, 'im')
                cur_ev_re = ev.value.re;
                cur_ev_im = ev.value.im;
            else
                cur_ev_re = NaN;
                cur_ev_im = NaN;
            end
            
            if isfield(ev, 'value_instab_2d') && any(ev.value_instab_2d < 0)
                unstable2d_re(end+1) = cur_ev_re;
                unstable2d_im(end+1) = cur_ev_im;
                reynolds_2d(end+1) = data.reynolds;
            end
            if isfield(ev, 'value_instab_3d') && any(ev.value_instab_3d < 0)
                unstable3d_re(end+1) = cur_ev_re;
                unstable3d_im(end+1) = cur_ev_im;
                reynolds_3d(end+1) = data.reynolds;
            end
        end
    end
    
    figure;
    hold on; grid on;
    
    if ~isempty(eigen_re)
        plot(eigen_re, eigen_im, 'kd', 'MarkerSize', 4, 'DisplayName', sprintf('Eigenvalues for Re %g, omega = %g+%gi', myRe, real(myOmega), imag(myOmega)));
    end
    
    plot(implied_alpha_2d_r, implied_alpha_2d_i, 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Implied alpha');
    
    all_re = unique([reynolds_2d, reynolds_3d]);
    cmap = colormap(jet(numel(all_re)));
    
    for k = 1:length(unstable2d_re)
        color_idx = find(all_re == reynolds_2d(k), 1);
        plot(unstable2d_re(k), unstable2d_im(k), 'bd', 'MarkerSize', 6, 'Color', cmap(color_idx,:),...
            'DisplayName', sprintf('Re = %g, 2D', reynolds_2d(k)));
    end
    
    for k = 1:length(unstable3d_re)
        color_idx = find(all_re == reynolds_3d(k), 1);
        plot(unstable3d_re(k), unstable3d_im(k), 'md', 'MarkerSize', 8, 'Color', cmap(color_idx,:),...
            'DisplayName', sprintf('Re = %g, 3D', reynolds_3d(k)));
    end

    xlabel('\alpha_{real}');
    ylabel('\alpha_{imag}');
    title(sprintf('Combined Plot: Specific Eigenvalues (Re = %g, omega = %g+%gi) and Unstable Modes\n(for omega_{real} in [%g, %g])', myRe, real(myOmega), imag(myOmega), minOmegaRe, maxOmegaRe));
    colormap(cmap);
    colorbar('Ticks', linspace(0, 1, numel(all_re)), 'TickLabels', all_re);
    hold off;
end
