%  http://www.k-wave.org/forum/topic/3d-simulations-of-lamb-wave-in-plate
clearvars;
% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================
scale = 1;
% define the source properties
source_freq = 10e6; % [Hz]
T_period = 1/(source_freq);
source_strength = 1e6; % [Pa]
source_cycles = 1; % number of tone burst cycles
source_position_y = 12;
% create the computational grid
PML_size = 10; % size of the PML in grid points
Nx = 32*scale - 2*PML_size; % number of grid points in the x direction
Ny = 1024*scale - 2*PML_size; % number of grid points in the y direction
% Nz = 300*scale - 2*PML_size; % number of grid points in the y direction
ppw = 5;
%rayleigh_speed = 2940;   % for Aluminium
rayleigh_speed = 3000;  % For steel
dx = rayleigh_speed/ (ppw * source_freq);   % [m]
%dx = 0.5e-3/scale; % grid point spacing in the x direction [m]
%dy = 0.5e-3/scale; % grid point spacing in the y direction [m]
dy=dx;
%kgrid = makeGrid(Nx, dx, Ny, dy);
kgrid = kWaveGrid(Nx, dx, Ny, dy);
% kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);
% define the medium properties
cp1 = 600; % compressional wave speed [m/s]
cs1 = 0; % shear wave speed [m/s]
rho1 = 200; % density [kg/m^3]
%alpha0_p1 = 5.0; % compressional wave absorption [dB/(MHz^2 cm)]
%alpha0_s1 = 5.0; % shear wave absorption [dB/(MHz^2 cm)]

%% medium properties for Steel
%
cp2 = 5960; % compressional wave speed [m/s]
cs2 = 3240; % shear wave speed [m/s]
rho2 = 7900; % density [kg/m^3]
%}

%% Medium properties for Aluminium
%{
cp2 = 6400; % compressional wave speed [m/s]
cs2 = 3150; % shear wave speed [m/s]
rho2 = 2700; % density [kg/m^3]
%}

%{
cp2 = 2656; % compressional wave speed [m/s]
cs2 = 1078; % shear wave speed [m/s]
rho2 = 1140; % density [kg/m^3]
%alpha0_p2 = 0.5; % compressional wave absorption [dB/(MHz^2 cm)]
%alpha0_s2 = 0.5; % shear wave absorption [dB/(MHz^2 cm)]
%}

% create the time array
cfl = 0.1;
%t_end = 800e-6;
t_end = (Nx*dx)*100.5/rayleigh_speed;
kgrid.t_array= makeTime(kgrid, cp2, cfl, t_end);
fs=1/kgrid.dt;

% define position of heterogeneous slab
% slab = ones(Nx, Ny, Nz);
% slab (1,:,:)=0 ;
% slab (19,:,:)=0 ;
slab = ones(Nx, Ny);
slab_strt = 4;
slab_end = Nx-3;
slab (1:slab_strt,:)=0 ;
slab (slab_end:end,:)=0 ;


%3D
% source_mask=zeros(Nx,Ny,Nz);
% source_mask(2,30,6)=1;
% source_mask(18,30,6)=1;
%2D
source_mask=zeros(Nx,Ny);
%ource_mask(slab_strt,source_position_y:source_position_y+5)=1;
source_mask(slab_strt,source_position_y)=1;

% assign the source
source.s_mask = source_mask;
%{
 source.sxx(1,:) = -source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);
 source.sxx(2,:) = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles);
%}

%
source.sxx(1,:) = source_strength*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%source.sxx(2,:) = source_strength*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%}

% define the sensor mask

sensor.mask = zeros(Nx, Ny);
%
sensors =  round(0.4*Ny):1:round(0.9*Ny);
source_sensor_dist = zeros(1,length(sensors));
i = 0; 
for k =  sensors
%for k = slab_strt:1:slab_end
%    sensor.mask(slab_strt,k:k+10) = 1;               %2D
    sensor.mask(slab_strt,k) = 1;               %2D
    %sensor.mask(k,round(0.8*Ny)) =1;
    i = i+1;
    source_sensor_dist(i) = (k - source_position_y)*dx*1e3;
end
%}
%sensor.mask(slab_strt, 0.8*Ny) = 1;

sensor.record = {'p','u'};
% set the input arguments
input_args = {'PMLSize', PML_size, 'PMLAlpha', 2, 'PlotPML', false, ...
'PMLInside', false, 'PlotScale', [-1, 1]*source_strength/2, ...
'DisplayMask', source.s_mask+sensor.mask  , 'DataCast',  'gpuArray-single'};

% =========================================================================
% ELASTIC SIMULATION
% =========================================================================
% define the medium properties
medium.sound_speed_compression = cp1*ones(Nx, Ny);
medium.sound_speed_compression(slab == 1) = cp2;
medium.sound_speed_shear = cs1*ones(Nx, Ny);
medium.sound_speed_shear(slab == 1) = cs2;
medium.density = rho1*ones(Nx, Ny);
medium.density(slab == 1) = rho2;
%{
medium.alpha_coeff_compression = alpha0_p1*ones(Nx, Ny);
medium.alpha_coeff_compression(slab == 1) = alpha0_p2;
medium.alpha_coeff_shear = alpha0_s1*ones(Nx, Ny);
medium.alpha_coeff_shear(slab == 1) = alpha0_s2;
%}

% run the elastic simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
save('lamb_2D.mat');

%% Plotting
% Scale the time and pressure values in SI unit
%
[t_sc, t_scale, t_prefix] = scaleSI(t_end);
[vec_sc, vec_scale, vec_prefix] = scaleSI(kgrid.x_vec);

p = sensor_data.p;
subplot(2,2,1)
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, p);
title('Preassure')
ux = sensor_data.ux;
subplot(2,2,2)
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, ux);
title('Ux')
uy = sensor_data.uy;
subplot(2,2,3)
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, uy);
title('Uy')
%colormap(getColorMap);
%save('SIM_Results2D.mat','p','ux','uy','fs','source_freq');
%{
figure;
no_of_sensors_m = length(sensor_data.p(:,1));
for m = 1:no_of_sensors_m
    subplot(no_of_sensors_m,1,m);
    plot(kgrid.t_array * t_scale,sensor_data.p(m,:));
%    if m < no_of_sensors_m, set(gca,'XTick',[]), end
    b = source_sensor_dist(m);
    legendinfo = [' Source-Sensor distance = ' num2str(b) 'mm'];
    grid on;
    legend(legendinfo);
end
%set(gca,'xtickMode', 'auto')
%xticks(0:50:t_end*t_scale);
sgtitle('Pressure over Time') 
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [Pa]');
grid on;
%}
%% Plottint the input signal
figure;
plot(kgrid.t_array*t_scale,source.sxx);
xlabel(['Time [' t_prefix 's]']);
%plot(source.sxx);
grid on;
xticks(0:2:t_end*t_scale);
%}

%% Plotting imagesc
%{
figure;
p=sensor_data.p;
subplot(2,2,1),imagesc(p),title('Preassure');
ux=sensor_data.ux;
subplot(2,2,2),imagesc(ux),title('Ux');
uy=sensor_data.uy;
subplot(2,2,3),imagesc(uy),title('Uy');
%}
