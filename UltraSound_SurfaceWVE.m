% Simulation of elastic wave in 3D solid medium. 
clearvars
% =========================================================================
                            %% SIMULATION
% =========================================================================
%% grid parameters

% define the driving signal
source_freq         = 2e6;    % [Hz]
source_strength     = 10e6;      % [Pa]
source_cycles       = 3;        % number of tone burst cycles
T_period = 1/(source_freq);
ppw       = 7;       % number of points per wavelength , use 1.5-3 for reyleigh wave
cfl   = 0.1;  
%rayleigh_speed = 4600;
rayleigh_speed = 3000;
dx = rayleigh_speed/ (ppw * source_freq);   % [m]
dy = dx;   
Nx = 128;           % number of grid points in the x (row) direction
%Ny = 8192;
Ny = 1024;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Assign material layer
PML_size        = 20;                      % [grid points]
sub_grid_x_strt = round(0.4*Nx);
sub_grid_x_end = round(0.9*Nx);
sub_grid_y_strt = round(0.05*Ny);
sub_grid_y_end = round(0.95*Ny);
source_position_x = round(sub_grid_x_strt - 40);
source_position_y = round(0.15*Ny);
sensor_position_x = round(sub_grid_x_strt - 20);
sensor_position_y = round(0.8*Ny);


% define position of heterogeneous slab
slab                = zeros(Nx, Ny);
%slab(sub_grid_x_strt:sub_grid_x_end, sub_grid_y_strt:sub_grid_y_end)   = 1;
slab(sub_grid_x_strt:end, sub_grid_y_strt:sub_grid_y_end)   = 1;

%% Medium properties
%
cp1                 = 1500;   	% compressional wave speed [m/s]
cs1                 = 0;        % shear wave speed [m/s]
rho1                = 1000;     % density [kg/m^3]
alpha0_p1           = 0.002;      % compressional absorption [dB/(MHz^2 cm)]
alpha0_s1           = 0.002;      % shear absorption [dB/(MHz^2 cm)]

%
% define the medium properties for the bottom layer ----- Stell
cp2                 = 5960;     % compressional wave speed [m/s]
cs2                 = 3240;     % shear wave speed [m/s]
rho2                = 7900;     % density [kg/m^3]
alpha0_p2           = 0;        % compressional absorption [dB/(MHz^2 cm)]
alpha0_s2           = 0;        % shear absorption [dB/(MHz^2 cm)]
%{

% define the medium properties for the bottom layer ----- Aluminum
cp2                 = 6361;     % compressional wave speed [m/s]
cs2                 = 3155;     % shear wave speed [m/s]
rho2                = 2709;     % density [kg/m^3]
alpha0_p2           = 0;        % compressional absorption [dB/(MHz^2 cm)]
alpha0_s2           = 0;        % shear absorption [dB/(MHz^2 cm)]
%

% define the medium properties for the bottom layer ----- Glass
cp2                 = 5700;     % compressional wave speed [m/s]
cs2                 = 3400;     % shear wave speed [m/s]
rho2                = 2600;     % density [kg/m^3]
alpha0_p2           = 0;        % compressional absorption [dB/(MHz^2 cm)]
alpha0_s2           = 0;        % shear absorption [dB/(MHz^2 cm)]
%}

%% Medium velocities

% define the medium properties
clear medium
medium.sound_speed_compression            = cp1*ones(Nx, Ny);
medium.sound_speed_compression(slab == 1) = cp2;
medium.sound_speed_shear                  = cs1*ones(Nx, Ny);
medium.sound_speed_shear(slab == 1)       = cs2;
medium.density                            = rho1*ones(Nx, Ny);
medium.density(slab == 1)                 = rho2;
medium.alpha_coeff_compression            = alpha0_p1*ones(Nx, Ny);
medium.alpha_coeff_compression(slab == 1) = alpha0_p2;
medium.alpha_coeff_shear                  = alpha0_s1*ones(Nx, Ny);
medium.alpha_coeff_shear(slab == 1)       = alpha0_s2;
%medium.alpha_power = 2;
%% Time Scale

t_end = 100e-6;
kgrid.makeTime(max(medium.sound_speed_compression(:),medium.sound_speed_shear(:)), cfl, t_end);


%% Sensor Mask

startpoint1 = [source_position_x,source_position_y];        % Source start position
%endpoint1 = [30,40];
startpoint2 = [sensor_position_x,sensor_position_y];    % sensor start position
%endpoint2 = [50,30];
angle1 = -pi*30/180;
angle2 = pi*30/180;
length_source = 30;
length_sensor = 3;
% generate the source geometry
source_mask = makeLine( Nx, Ny, startpoint1, angle1, length_source);

% assign the source
source.s_mask = source_mask;
source.sxx = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);
source.syy = source.sxx;

sensor.mask = zeros(Nx,Ny);
%sensor.mask = makeLine( Nx, Ny, startpoint2, angle2, lengt_h2);
%sensor.mask(sub_grid_strt - 10, sensor_position_y) = 1;

%
%sensors =  round(31e-3/dx):round(2e-3/dx):round(33e-3/dx);
%sensors = round(linspace(round(0.3*Ny),round(0.9*Ny),12));
sensors = round(0.4*Ny):round(10e-3/dx):round(0.4*Ny)+9*round(10e-3/dx);
source_sensor_dist = zeros(1,length(sensors));
i = 0; 
for k = sensors
%for k = slab_strt:1:slab_end
%    sensor.mask(slab_strt,k:k+10) = 1;               %2D
%    sensor.mask(sub_grid_x_strt,k) = 1;               %2D
    sensor.mask(sensor_position_x , k) = 1;
    %sensor.mask(k,round(0.8*Ny)) =1;
    i = i+1;
    source_sensor_dist(i) = (k - (source_position_y - round(length_source/2)))*dx;
end
%
%}

%% define a custom display mask showing the position of the interface from

display_mask = false(Nx, Ny);
%display_mask(Nx/2, :) = 1;
%display_mask(slab==1) = 1;
display_mask(sub_grid_x_strt, :) = 1;

%% Data Record
%sensor.record = {'u','u_split_field'};
sensor.record = {'p','u'};
% define input arguments
input_args = {'PlotScale',[-1, 1]*source_strength,'PlotPML', true,'PMLSize', PML_size,  'DisplayMask',  display_mask + source.s_mask + sensor.mask,'DataCast', 'gpuArray-single', 'PMLAlpha',4};
%{
input_args = {  'PlotScale',[-1, 1]*source_strength,...
                'PlotPML', false,...
                'RecordMovie', true,...
                'DisplayMask', source.s_mask+sensor.mask,...
                'MovieName','surface_2d_scale5_pint_sc_1scale',...
                'MovieArgs', {'FrameRate', 10},...
                'DataCast', 'gpuArray-single',...
                'PlotFreq', 200,...
                'MovieProfile', 'MPEG-4', ...
                'PMLAlpha',4    }; 
%}

%% run the simulation

sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
save('US_stl.mat');

% =========================================================================
                               %% VISUALISATION
% =========================================================================
%
%[t_sc, t_scale, t_prefix] = scaleSI(t_end);
[p_sc, p_scale, p_prefix] = scaleSI(sensor_data.p(1,:));
[d_sc, d_scale, d_prefix] = scaleSI( source_sensor_dist(end));
[t_sc, t_scale, t_prefix] = scaleSI(t_end);

%% For plotting multiple sensor point time series data
%
figure;
no_of_sensors_m = length(sensor_data.p(:,1));
for m = 1:no_of_sensors_m
    subplot(no_of_sensors_m,1,m);
   plot(kgrid.t_array * t_scale,(sensor_data.p(m,:))*p_scale);
%    plot(kgrid.t_array * t_scale,(sensor_data.p(m,:)/max(sensor_data.p(m,:))));
%    if m < no_of_sensors_m, set(gca,'XTick',[]), end
    b = source_sensor_dist(m)*d_scale;
    grid on;
    legend( [' Source-Sensor distance = ' num2str(b)  d_prefix 'm'])
end
%set(gca,'xtickMode', 'auto')
%xticks(0:50:t_end*t_scale);
title('Pressure over Time') 
xlabel(['Time [' t_prefix 's]']);
ylabel(['Pressure [' p_prefix ' Pa]']);
%grid on;
%}

figure;
q = waterfall(kgrid.t_array * t_scale,1:10,(sensor_data.p)*p_scale);         % plotting the data in waterfall
q.EdgeColor = 'b';
title('Pressure over Time') 
xlabel('Time[µs]');
ylabel('Position[mm]');
zlabel('Amplitude [KPa]');

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
%}