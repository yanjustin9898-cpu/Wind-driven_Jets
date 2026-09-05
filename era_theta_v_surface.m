% Script: calc_era5_surface_theta_v.m
% Calculates surface theta_v from ERA5 single level NetCDF dataset

filename = 'ATLera5_monthly_temp_2005-2009_surface.nc';

% --- 1. Read Variables ---
try
    t2m = double(ncread(filename, 't2m')); % 2m temperature [K]
    d2m = double(ncread(filename, 'd2m')); % 2m dewpoint temperature [K]
    sp  = double(ncread(filename, 'sp'));  % Surface pressure [Pa]
catch
    t2m = double(ncread(filename, '2m_temperature'));
    d2m = double(ncread(filename, '2m_dewpoint_temperature'));
    sp  = double(ncread(filename, 'surface_pressure'));
end

% --- 2. Constants ---
P0 = 1000.0;    % Reference pressure [hPa]
kappa = 0.2854; % Rd / cp for dry air

% --- 3. Compute Vapor Pressure e [hPa] ---
T_d_C = d2m - 273.15; % Convert K to C
e = 6.112 .* exp((17.67 .* T_d_C) ./ (T_d_C + 243.5));

% --- 4. Convert Surface Pressure to hPa ---
P_hpa = sp ./ 100.0;

% --- 5. Compute Mixing Ratio r [kg/kg] ---
r = 0.622 .* (e ./ (P_hpa - e));

% --- 6. Compute Surface Virtual Temperature Tv [K] ---
Tv_sfc = t2m .* (1.0 + 0.61 .* r);

% --- 7. Compute Surface Virtual Potential Temperature theta_v [K] ---
theta_v_sfc = Tv_sfc .* ((P0 ./ P_hpa) .^ kappa);

fprintf('Successfully computed surface theta_v.\n');
fprintf('Array dimensions: %s\n', mat2str(size(theta_v_sfc)));

save('ATLera5_theta_v_surface2m.mat', 'theta_v_sfc', '-v7.3');
