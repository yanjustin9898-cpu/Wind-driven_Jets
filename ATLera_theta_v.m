% Script: calculate_era5_theta_v.m
% Calculates Virtual Potential Temperature (theta_v) from ERA5 NetCDF data

filename = 'ATLera5_monthly_temp_2005-2009_1000-700hPa.nc';

% --- 1. Read Variables from NetCDF ---
% Note: Standard ERA5 NetCDF names are 't' (temp), 'q' (humidity), 'level' (pressure)
try
    T = ncread(filename, 't');          % Temperature [K]
    q = ncread(filename, 'q');          % Specific humidity [kg/kg]
    level = ncread(filename, 'pressure_level');  % Pressure levels [hPa]
catch
    % Alternate variable names depending on download source (e.g., CDS beta)
    T = ncread(filename, 'temperature');
    q = ncread(filename, 'specific_humidity');
    level = ncread(filename, 'isobaricInhPa');
end

% Ensure variables are double precision for math operations
T = double(T);
q = double(q);
p = double(level);

% --- 2. Constants ---
P0 = 1000.0;    % Reference pressure [hPa]
kappa = 0.2854; % Dry air gas constant / specific heat (Rd / cp)

% --- 3. Compute Water Vapor Mixing Ratio (r) ---
r = q ./ (1.0 - q); % [kg/kg]

% --- 4. Compute Virtual Temperature (Tv) ---
Tv = T .* (1.0 + 0.61 .* r); % [K]

% --- 5. Compute Virtual Potential Temperature (theta_v) ---
% ERA5 array structure is typically [longitude, latitude, level, time].
% Reshape 'level' into a 3D/4D vector [1 x 1 x nlevels x 1] for implicit expansion.
p_grid = reshape(p, 1, 1, []);

% Apply Poisson's equation
theta_v = Tv .* ((P0 ./ p_grid) .^ kappa); % [K]

% --- 6. Save or Export Results ---
% Display confirmation and dataset info
fprintf('Successfully calculated theta_v.\n');
fprintf('Array dimensions: %s\n', mat2str(size(theta_v)));

geo   = double(ncread(filename, 'z'));     % Geopotential [m^2/s^2]
lon = ncread(filename, 'longitude');
lat = ncread(filename, 'latitude');
time = ncread('era5_monthly_temp_2005-2009_1000-700hPa.nc', 'valid_time');

% --- 3. Convert Geopotential to Geopotential Height Z [m] ---
g0 = 9.80665;    % Standard gravity [m/s^2]
Z = geo ./ g0;   % Height in meters above sea level [lon, lat, level, time]

% --- 4. Define Target Regular Height Grid ---
% 1000 hPa to 700 hPa spans roughly 0 to 3000 m ASL
z_target = 0:100:3000; % [m] vertical grid from 0m to 3000m at 100m steps

% --- 5. Vertical Interpolation (Pressure -> Height Coordinates) ---
[nlon, nlat, nlevel, ntime] = size(theta_v);
n_target = length(z_target);

% Reshape 4D arrays to 2D [nlevel x N] where columns are individual vertical profiles
Z_profiles  = reshape(permute(Z, [3, 1, 2, 4]), nlevel, []);
th_profiles = reshape(permute(theta_v, [3, 1, 2, 4]), nlevel, []);

N = size(Z_profiles, 2);
th_regrid_col = NaN(n_target, N);

% Interpolate each column independently
for i = 1:N
    z_col  = Z_profiles(:, i);
    th_col = th_profiles(:, i);

    % Sort profile monotonically by height (required by interp1)
    [z_sorted, idx] = sort(z_col);
    th_sorted = th_col(idx);

    % Remove NaNs if present
    valid = ~isnan(z_sorted) & ~isnan(th_sorted);

    if sum(valid) >= 2
        % Linear 1D interpolation; points outside the column range return NaN
        th_regrid_col(:, i) = interp1(z_sorted(valid), th_sorted(valid), z_target, 'linear', NaN);
    end
end

% --- 6. Reshape Output Array Back to 4D [nlon x nlat x nheight x ntime] ---
theta_v_height = permute(reshape(th_regrid_col, [n_target, nlon, nlat, ntime]), [2, 3, 1, 4]);

ERA5_ThetaV.theta_v   = theta_v_height; % 4D matrix [lon x lat x height x time]
ERA5_ThetaV.longitude = lon;            % 1D vector
ERA5_ThetaV.latitude  = lat;            % 1D vector
ERA5_ThetaV.height    = z_target;       % 1D vector (meters ASL)
ERA5_ThetaV.time      = time;           % 1D vector

% Save with -v7.3 for large 4D datasets (> 2 GB support)
save('ATLera5_theta_v_height_coords.mat', 'ERA5_ThetaV', '-v7.3');

