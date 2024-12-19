% ------------------------------ %
%    Action Potential Simulation %
% ------------------------------ %
% This script simulates neuronal action potentials (AP) using the
% Hodgkin-Huxley model augmented with BK (Big Potassium) channels.
% It explores how varying BK conductance (gBK) and calcium concentration (Ca)
% affect the AP waveform, particularly the half-width of the AP.
% Additionally, it visualizes the half-width directly on the AP traces.

close all; clear; clc;

%% ---------------------- %%
%%   Simulation Setup     %%
%% ---------------------- %%

% Time parameters
dt = 0.01;          % Time step (ms)
tmax = 50;          % Maximum time (ms)
t = 0:dt:tmax;      % Time vector

% Reversal potentials (mV)
ENa = 115;   % Sodium reversal potential
EK = -12;    % Potassium reversal potential
El = 10.6;   % Leak reversal potential
EBK = -80;   % BK channel reversal potential

% Maximum conductances (mS/cm^2)
gNa = 120;   % Sodium conductance
gK = 36;     % Potassium conductance
gL = 0.3;    % Leak conductance

% Capacitance (uF/cm^2)
Cm = 1.0;

% BK channel parameters
gBK_values = [1, 10, 20];      % Different BK conductances (mS/cm^2)
Ca_values = [0.5];              % Calcium concentrations (nM)

% Initialize storage for membrane potentials and half-widths
V_conditions = zeros(length(t), numel(gBK_values), numel(Ca_values));
half_widths = zeros(numel(gBK_values), numel(Ca_values));
rise_times = zeros(numel(gBK_values), numel(Ca_values));
fall_times = zeros(numel(gBK_values), numel(Ca_values));

%% ------------------------ %%
%%   Rate Function Setup    %%
%% ------------------------ %%

% Define a small epsilon to handle singularities
epsilon = 1e-6;

% Anonymous functions to compute rate constants with singularity handling

% Function to compute alpha_m with handling for V approaching 25 mV
compute_alpha_m = @(V) ...
    (abs(25 - V) < epsilon) * 1 + ...  % Limit as V approaches 25
    (0.1*(25 - V)./(exp((25 - V)/10) - 1)) .* (abs(25 - V) >= epsilon);

% Function to compute alpha_n with handling for V approaching 10 mV
compute_alpha_n = @(V) ...
    (abs(10 - V) < epsilon) * 0.1 + ...  % Limit as V approaches 10
    (0.01*(10 - V)./(exp((10 - V)/10) - 1)) .* (abs(10 - V) >= epsilon);

% Function to compute alpha_b with handling for V approaching 25 mV
compute_alpha_b = @(V, Ca_dep) ...
    0.02 * ( (abs(25 - V) < epsilon) * 1 + ...  % Limit as V approaches 25
    ((25 - V)./(exp((25 - V)/10) - 1)) ) .* (abs(25 - V) >= epsilon) .* Ca_dep;

%% ------------------------- %%
%%      Simulation Loop      %%
%% ------------------------- %%

% Loop through different BK conductances
for c_gBK = 1:numel(gBK_values)
    for c_Ca = 1:numel(Ca_values)
        % Current condition parameters
        gBK = gBK_values(c_gBK);
        Ca = Ca_values(c_Ca);
        
        % Calculate calcium-dependent factor
        Ca_dep = Ca / (Ca + 100);
        
        % Initialize variables
        V = -70 * ones(size(t));   % Membrane potential (mV)
        
        % Initialize gating variables to steady-state at V = -70 mV
        % Compute alpha and beta at V = -70 mV
        alpha_m_inf = compute_alpha_m(V(1));
        beta_m_inf = 4 * exp(-V(1)/18);
        m_inf = alpha_m_inf / (alpha_m_inf + beta_m_inf);
        
        alpha_h_inf = 0.07 * exp(-V(1)/20);
        beta_h_inf = 1 / (exp((30 - V(1))/10) + 1);
        h_inf = alpha_h_inf / (alpha_h_inf + beta_h_inf);
        
        alpha_n_inf = compute_alpha_n(V(1));
        beta_n_inf = 0.125 * exp(-V(1)/80);
        n_inf = alpha_n_inf / (alpha_n_inf + beta_n_inf);
        
        alpha_b_inf = compute_alpha_b(V(1), Ca_dep);
        beta_b_inf = 0.1 * exp(-V(1)/20);
        b_inf = alpha_b_inf / (alpha_b_inf + beta_b_inf);
        
        % Initialize gating variables to steady-state values
        m = m_inf * ones(size(t));
        h = h_inf * ones(size(t));
        n = n_inf * ones(size(t));
        b = b_inf * ones(size(t));
        
        % Define input current
        I = zeros(size(t));
        I(t > 10 & t < 11) = 20;  % Current injection between 10 ms and 11 ms (uA/cm^2)
        
        % Simulation using Euler's method
        for i = 1:length(t)-1
            % Compute rate constants
            alpha_m = compute_alpha_m(V(i));
            beta_m = 4 * exp(-V(i)/18);
            
            alpha_h = 0.07 * exp(-V(i)/20);
            beta_h = 1 / (exp((30 - V(i))/10) + 1);
            
            alpha_n = compute_alpha_n(V(i));
            beta_n = 0.125 * exp(-V(i)/80);
            
            alpha_b = compute_alpha_b(V(i), Ca_dep);
            beta_b = 0.1 * exp(-V(i)/20);
            
            % Update gating variables
            m(i+1) = m(i) + dt * (alpha_m * (1 - m(i)) - beta_m * m(i));
            h(i+1) = h(i) + dt * (alpha_h * (1 - h(i)) - beta_h * h(i));
            n(i+1) = n(i) + dt * (alpha_n * (1 - n(i)) - beta_n * n(i));
            b(i+1) = b(i) + dt * (alpha_b * (1 - b(i)) - beta_b * b(i));
            
            % Calculate ionic currents
            INa = gNa * m(i)^3 * h(i) * (V(i) - ENa);
            IK = gK * n(i)^4 * (V(i) - EK);
            IL = gL * (V(i) - El);
            IBK = gBK * b(i) * (V(i) - EBK);
            
            % Total ionic current
            I_ion = INa + IK + IL + IBK;
            
            % Update membrane potential
            V(i+1) = V(i) + dt * (-I_ion + I(i)) / Cm;
        end

        % Store membrane potential for current condition
        V_conditions(:, c_gBK, c_Ca) = V;
        
        % Calculate and store half-width and crossing times
        [half_widths(c_gBK, c_Ca), rise_times(c_gBK, c_Ca), fall_times(c_gBK, c_Ca)] = calculate_half_width(t, V);
    end
end

%% --------------------- %%
%%        Plotting        %%
%% --------------------- %%

% Plotting Action Potentials with Half-Width Annotations
figure('Name', 'Action Potentials with Half-Width Annotations');
hold on;
colors = {'r', 'g', 'b'};
styles = {'-', '--', ':'};

for c_gBK = 1:numel(gBK_values)
    for c_Ca = 1:numel(Ca_values)
        V = V_conditions(:, c_gBK, c_Ca);
        [half_width, rise_time, fall_time] = deal(half_widths(c_gBK, c_Ca), rise_times(c_gBK, c_Ca), fall_times(c_gBK, c_Ca));
        
        % Plot AP trace
        plot(t, V, 'Color', colors{c_gBK}, 'LineStyle', styles{c_Ca}, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('gBK=%d mS/cm², Ca=%.1f nM', gBK_values(c_gBK), Ca_values(c_Ca)));
         
        % Calculate half-maximum voltage
        V_peak = max(V);
        half_V = V_peak / 2;
        
        % Plot horizontal line at half_V
        plot([min(t) max(t)], [half_V half_V], 'k--', 'LineWidth', 0.5);
        
        % Plot vertical lines at rise_time and fall_time
        plot([rise_time rise_time], [min(V) half_V], 'k:', 'LineWidth', 1);
        plot([fall_time fall_time], [min(V) half_V], 'k:', 'LineWidth', 1);
        
        % Plot markers at rise and fall points
        plot(rise_time, half_V, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        plot(fall_time, half_V, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        
        % Annotate half-width
        text((rise_time + fall_time)/2, half_V + 5, sprintf('Half-Width: %.2f ms', half_width), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', colors{c_gBK});
    end
end

% Customize AP Plot
legend('Location','best');
title('Action Potentials with Half-Width Annotations');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;
hold off;

% Plotting Half-Width vs BK Conductance
figure('Name', 'AP Half-Width vs BK Conductance');
hold on;
markers = {'o-', 's-', 'd-'};
for c_gBK = 1:numel(gBK_values)
    plot(gBK_values(c_gBK), half_widths(c_gBK,1), markers{c_gBK}, 'LineWidth', 2, 'Color', colors{c_gBK}, ...
         'DisplayName', sprintf('gBK=%d mS/cm², Ca=%.1f nM', gBK_values(c_gBK), Ca_values(1)));
end
xlabel('BK Conductance (mS/cm²)');
ylabel('Action Potential Half-Width (ms)');
title(['AP Half-Width vs BK Conductance (Ca = ', num2str(Ca_values(1)), ' nM)']);
grid on;
legend('Location','best');
hold off;

%% -------------------------------------- %%
%%    Half-Width Calculation Function    %%
%% -------------------------------------- %%

function [half_width, rise_time, fall_time] = calculate_half_width(t, V)
    %CALCULATE_HALF_WIDTH Calculates the half-width of an action potential.
    %
    % Parameters:
    %   t - Time vector (ms)
    %   V - Membrane potential vector (mV)
    %
    % Returns:
    %   half_width - Half-width of the action potential (ms)
    %   rise_time  - Time when V crosses half_V on rising phase (ms)
    %   fall_time  - Time when V crosses half_V on falling phase (ms)

    % Find peak voltage and its index
    [V_peak, peak_idx] = max(V);
    half_V = V_peak / 2;
    
    % Find the time point where V crosses half_V on the rising phase
    rising_phase = V(1:peak_idx);
    rising_indices = find(rising_phase >= half_V, 1, 'first');
    if isempty(rising_indices)
        rise_time = t(peak_idx);
    else
        % Linear interpolation for more accurate crossing time
        if rising_indices == 1
            rise_time = t(rising_indices);
        else
            V1 = V(rising_indices - 1);
            V2 = V(rising_indices);
            t1 = t(rising_indices - 1);
            t2 = t(rising_indices);
            rise_time = t1 + (half_V - V1) * (t2 - t1) / (V2 - V1);
        end
    end
    
    % Find the time point where V crosses half_V on the falling phase
    falling_phase = V(peak_idx:end);
    falling_indices = find(falling_phase <= half_V, 1, 'first');
    if isempty(falling_indices)
        fall_time = t(end);
    else
        idx = peak_idx + falling_indices - 1;
        if falling_indices == 1
            fall_time = t(idx);
        else
            V1 = V(idx - 1);
            V2 = V(idx);
            t1 = t(idx - 1);
            t2 = t(idx);
            fall_time = t1 + (half_V - V1) * (t2 - t1) / (V2 - V1);
        end
    end
    
    % Calculate half-width
    half_width = fall_time - rise_time;
end
