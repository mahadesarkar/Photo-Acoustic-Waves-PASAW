% Simulation of elastic wave in 3D solid medium. 
clearvars
% =========================================================================
                            %% SIMULATION
% =========================================================================
%% grid parameters
source_f0 = 10e6; % [Hz]
T_period = 1/(source_f0);
ppw       = 7;                   % number of points per wavelength , use 1.5-3 for reyleigh wave
cfl   = 0.1;  
%rayleigh_speed = 4600;
rayleigh_speed = 3000;
dx = rayleigh_speed/ (ppw * source_f0);   % [m]
dy = dx;   
dz = dx;
%Nx = 128;                       % number of grid points in the x (row) direction
Nx = 256;
Ny = 4096;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Assign material layer
sub_grid_strt = round(0.4*Nx);
source_position_x = sub_grid_strt;
source_position_y = 40;

% Assign crack dimension
%{
crack_start_x = sub_grid_strt;
crack_end_x = sub_grid_strt + 20;
crack_start_y = round(0.60*Ny);
crack_end_y = round(0.61*Ny);
%}

%% Medium properties
%
% Velocity & density for Steel medium layer
medium.sound_speed_compression = 800 * ones(Nx, Ny);   % 9138[m/s]
medium.sound_speed_shear       = 500 * ones(Nx, Ny);   % 4675 [m/s]
medium.density                 = 200 * ones(Nx, Ny);   % [kg/m^3]
medium.alpha_coeff_compression           = 0.002* ones(Nx, Ny);      % compressional absorption [dB/(MHz^2 cm)]
medium.alpha_coeff_shear         = 0.002* ones(Nx, Ny);      % shear absorption [dB/(MHz^2 cm)]
%}

%{
% Vlocity & density for Si medium layer
medium.sound_speed_compression = 930 * ones(Nx, Ny);   % 9138[m/s] 
medium.sound_speed_shear       = 540 * ones(Nx, Ny);   % 4675 [m/s]
medium.density                 = 240 * ones(Nx, Ny);   % [kg/m^3]
%}

% define the properties of the upper layer of the propagation medium
% For steel medium
%
medium.sound_speed_compression(sub_grid_strt:end,:) = 5960;   % [m/s]
medium.sound_speed_shear(sub_grid_strt:end,:)       = 3240;   % [m/s]
medium.density(sub_grid_strt:end,:)                 = 7900;   % [kg/m^3]
medium.alpha_coeff_compression(sub_grid_strt:end,:)           = 0;        % compressional absorption [dB/(MHz^2 cm)]
medium.alpha_coeff_shear(sub_grid_strt:end,:)           = 0;        % shear absorption [dB/(MHz^2 cm)]
%}
%{
% For Aluminium medium
medium.sound_speed_compression(sub_grid_strt:end,:) = 6400;   % [m/s]
medium.sound_speed_shear(sub_grid_strt:end,:)       = 3150;   % [m/s]
medium.density(sub_grid_strt:end,:)                 = 2700;   % [kg/m^3]
%}
%{
% For Si medium
medium.sound_speed_compression(sub_grid_strt:end,:) = 9300;   % [m/s]
medium.sound_speed_shear(sub_grid_strt:end,:)       = 5800;   % [m/s]
medium.density(sub_grid_strt:end,:)                 = 2400;   % [kg/m^3]
%}

% Crack medium properties and dimension defination
%{
medium.sound_speed_compression(crack_start_x:crack_end_x,crack_start_y:crack_end_y) = 800;   % [m/s]
medium.sound_speed_shear(crack_start_x:crack_end_x,crack_start_y:crack_end_y)       = 500;   % [m/s]
medium.density(crack_start_x:crack_end_x,crack_start_y:crack_end_y)                 = 200;   % [kg/m^3]
%}

% define the absorption properties
%medium.alpha_coeff_compression = 0.1 ; %0.1; % [dB/(MHz^2 cm)]
%medium.alpha_coeff_shear       = 0.5;   %0.5; % [dB/(MHz^2 cm)]
%medium.alpha_power = 2;      % Not erquired for Pstd_3d medium defination.
%medium.BonA =0.5;

%% Time Scale

%Nt = 1000;
%dt = 1e-5;
%kgrid.setTime(Nt, dt);
%t_end = (Nx*dx)*5/rayleigh_speed;   % [s]
t_end = 100e-6;
kgrid.makeTime(max(medium.sound_speed_compression(:),medium.sound_speed_shear(:)), cfl, t_end);
time = kgrid.t_array(kgrid.t_array<1*T_period);              %  5 cycle time array for sinusoidal input pulse
[t_sc, t_scale, t_prefix] = scaleSI(t_end);

%% Pressure from temperature distribution
Pulse_duraton = 10e-9;
center = 10e-9;
laser_energy = 0.1e-3;
power = laser_energy/Pulse_duraton;
area = pi*(0.5e-3)^2;
intensity = power/area;


%  tem = temperature_time_gaussian_pulse(intensity_max,difusivity,conductivity,t0,tao,total_time)
% call tem.pressure(1,:)

% Diffusivity in m²/s and conductivity in W/m/K, 1cm-1 = 100 m-1, thermal
% conductivity from Scruby page 227

tem = temp_init_pres(intensity,1.3e-5,50,0.5e-3,center,Pulse_duraton);
%tem = temp_init_pres(1.6e10,0.88e-4,150e3,1e-3,20e-9,10e-9);
press = squeeze(tem.pressure);   
press_length = length(press(:,1));
source_mag = max(press(press_length/2,:));
%plot_time = kgrid.t_array(kgrid.t_array<= 30e-9);

%% Stress source DEfination

source.s_mask = zeros(Nx, Ny);
%source.s_mask= makeDisc(Nx,Ny, source_position_x, source_position_y, 3);
source.s_mask(source_position_x, source_position_y:source_position_y+(press_length-1)) = 1;
source.s_mode = 'dirichlet';

%max_pressure = max(tem.pressure(1,:));
%source.sxx = max_pressure*toneBurst(1/kgrid.dt, source_f0, 2);                %  Toneburst signal
%source.sxx = source_mag*exp(-((kgrid.t_array - T_period)/(T_period/4)).^2);   %  Gaussian signal
%source.sxx = 5*tem.pressure(1,:);  % Single initial pressure time series
%source.sxx = 0.5*max_pressure;
%source_mag = max(source.sxx);
%}

% Assign the initial pressure time series as pressure line source
%
for initialp_no = 1:press_length
    source.sxx(initialp_no,:) = press(initialp_no,2:end);     % Value scalled to 5 times for better resolution
end
%}
source.syy = source.sxx;
%source.sxy = 0;

%% Sensor Mask
%  Uncomment this line to apply the source for simulation

sensor.mask = zeros(Nx,Ny);
sensors = round(0.3*Ny):round(2e-3/dx):round(0.3*Ny)+9*round(2e-3/dx);
source_sensor_dist = zeros(1,length(sensors));
i = 0; 
for k = sensors
%for k = slab_strt:1:slab_end
%    sensor.mask(slab_strt,k:k+10) = 1;               %2D
%    sensor.mask(sub_grid_strt - 30,k) = 1;               %2D
    sensor.mask(sub_grid_strt -30,k) = 1;               %2D
    %sensor.mask(k,round(0.8*Ny)) =1;
    i = i+1;
    source_sensor_dist(i) = (k - (source_position_y + length(press_length)))*dx;
end
%

%% define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
%display_mask(sub_grid_strt, :) = 1;
%display_mask(crack_start_x:crack_end_x,crack_start_y:crack_end_y) = 1;         % Display the crack


%% Data Record
%sensor.record = {'u','u_split_field'};
sensor.record = {'p','u'};
% define input arguments
%
%input_args = {'PlotScale', [-1, 1] , 'PlotPML', false,...
%    'DisplayMask', display_mask,'MovieName','elastic_2d', 'DataCast', 'gpuArray-single','PMLAlpha',4};
input_args = {'PlotScale',[-source_mag, source_mag],'PlotPML', true, 'DisplayMask', display_mask + source.s_mask + sensor.mask, 'DataCast', 'single','PMLAlpha',4}; 
%{
input_args = {  'PlotScale', [-source_mag, source_mag],...
                'PlotPML', false,...
                'RecordMovie', true,...
                'DisplayMask', display_mask+source.s_mask+sensor.mask,...
                'MovieName','surface_2d_scale5_pint_sc_1scale',...
                'MovieArgs', {'FrameRate', 10},...
                'PlotFreq', 200,...
                'MovieProfile', 'MPEG-4', ...
                'DataCast', 'gpuArray-single',...
                'PMLAlpha',4    }; 
%}
%% run the simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
save('Laser_stl.mat');
%}
% =========================================================================
                               %% VISUALISATION
% =========================================================================
%
%[t_sc, t_scale, t_prefix] = scaleSI(t_end);
[p_sc, p_scale, p_prefix] = scaleSI(sensor_data.p(1,:));
[d_sc, d_scale, d_prefix] = scaleSI( source_sensor_dist(end));
    
%% For plotting multiple sensor point time series data
%
figure;
no_of_sensors_m = length(sensor_data.p(:,1));
for m = 1:no_of_sensors_m
    subplot(no_of_sensors_m,1,m);
    plot(kgrid.t_array * t_scale,sensor_data.p(m,:)*p_scale);
%    if m < no_of_sensors_m, set(gca,'XTick',[]), end
    b = source_sensor_dist(m)*d_scale;
    grid on;
    legend( [' Source-Sensor distance = ' num2str(b)  d_prefix 'm'])
end
%set(gca,'xtickMode', 'auto')
%xticks(0:50:t_end*t_scale);
sgtitle('Pressure over Time') 
xlabel(['Time [' t_prefix 's]']);
ylabel(['Pressure [' p_prefix ' Pa]']);
%grid on;
%}
%% SOurce and sensor mask plot 
%
figure;
[vec_sc, vec_scale, vec_prefix] = scaleSI(kgrid.x_vec);
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, ...
    double(source.s_mask | sensor.mask | display_mask), [-1 1]);
colormap(getColorMap);
ylabel(['x position [' vec_prefix 'm]']);
xlabel(['yposition [' vec_prefix 'm]']);
axis image;


figure;
q = waterfall(kgrid.t_array * t_scale,1:10,(sensor_data.p)*p_scale);         % plotting the data in waterfall
q.EdgeColor = 'b';
title('Pressure over Time') 
xlabel('Time[µs]');
ylabel('Position[mm]');
zlabel('Amplitude [KPa]');


%}