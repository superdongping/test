close all;
% Time
dt = 0.01;
tmax = 50;
t = 0:dt:tmax;

% Constants
ENa = 115;
EK = -12;
El = 10.6;
EBK = -80;

gNa = 120;
gK = 36;
gL = 0.3;

Cm = 1.0;

% Different BK conductance and Calcium levels
gBK_values = [1, 10, 20];  % Low, middle, high
Ca_values = [0.5]; % Low, medium, high Ca concentration in nM

% Initialize variables for storage
V_conditions = zeros(length(t), numel(gBK_values), numel(Ca_values));

% Loop through different conditions
for c_gBK = 1:numel(gBK_values)
    for c_Ca = 1:numel(Ca_values)
        gBK = gBK_values(c_gBK);
        Ca = Ca_values(c_Ca);
        
        % Initialize variables
        V = -70*ones(size(t));
        m = zeros(size(t));
        h = zeros(size(t));
        n = zeros(size(t));
        b = zeros(size(t));
        
        I = zeros(size(t));
        I(t > 10 & t < 11) = 20;

        % Simulation
        for i = 1:length(t)-1
            alpha_m = 0.1*(25-V(i))/(exp((25-V(i))/10)-1);
            beta_m = 4*exp(-V(i)/18);
            m(i+1) = m(i) + dt*(alpha_m*(1-m(i)) - beta_m*m(i));
            
            alpha_h = 0.07*exp(-V(i)/20);
            beta_h = 1/(exp((30-V(i))/10)+1);
            h(i+1) = h(i) + dt*(alpha_h*(1-h(i)) - beta_h*h(i));
            
            alpha_n = 0.01*(10-V(i))/(exp((10-V(i))/10)-1);
            beta_n = 0.125*exp(-V(i)/80);
            n(i+1) = n(i) + dt*(alpha_n*(1-n(i)) - beta_n*n(i));
            
            Ca_dep = Ca / (Ca + 100);
            alpha_b = 0.02*(25-V(i))/(exp((25-V(i))/10)-1) * Ca_dep;
            beta_b = 0.1*exp(-V(i)/20);
            b(i+1) = b(i) + dt*(alpha_b*(1-b(i)) - beta_b*b(i));
            
            INa = gNa * m(i)^3 * h(i) * (V(i) - ENa);
            IK = gK * n(i)^4 * (V(i) - EK);
            IL = gL * (V(i) - El);
            IBK = gBK * b(i) * (V(i) - EBK);
            
            I_ion = INa + IK + IL + IBK;
            V(i+1) = V(i) + dt*(-I_ion + I(i)) / Cm;
        end

        % Save this condition
        V_conditions(:, c_gBK, c_Ca) = V;
    end
end

% Plotting
figure;
hold on;
colors = {'r', 'g', 'b'};
styles = {'-', '--', ':'};

for c_gBK = 1:numel(gBK_values)
    for c_Ca = 1:numel(Ca_values)
        plot(t, V_conditions(:, c_gBK, c_Ca), ...
             'Color', colors{c_gBK}, 'LineStyle', styles{c_Ca}, ...
             'DisplayName', sprintf('gBK=%d, Ca=%d', gBK_values(c_gBK), Ca_values(c_Ca)));
    end
end

legend('Location','best');
title('Action Potentials with Different BK Conductances and Ca Concentrations');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
hold off;
