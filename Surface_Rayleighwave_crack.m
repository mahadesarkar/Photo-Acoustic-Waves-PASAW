% Simulation of elastic wave in 3D solid medium. 
clearvars
% =========================================================================
                            %% SIMULATION
% =========================================================================
%% grid parameters

source_f0 = 10e6; % [Hz]
T_period = 1/(source_f0);
ppw       = 30;       % number of points per wavelength , use 1.5-3 for reyleigh wave
cfl   = 0.1;  
%rayleigh_speed = 4600;
rayleigh_speed = 3000;
dx = rayleigh_speed/ (ppw * source_f0);   % [m]
dy = dx;   
dz = dx;
Nx = 128;           % number of grid points in the x (row) direction
%Ny = 8192;
Ny = 2048;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Assign material layer
sub_grid_strt = round(0.2*Nx);
source_position_x = sub_grid_strt;
source_position_y = 30;


% Assign crack dimension
crack_start_x = sub_grid_strt+5;
crack_end_x = sub_grid_strt + 20;
crack_start_y = round(0.60*Ny);
crack_end_y = round(0.61*Ny);

%% Medium properties
%
% Velocity & density for Steel medium layer
medium.sound_speed_compression = 800 * ones(Nx, Ny);   % 9138[m/s]
medium.sound_speed_shear       = 500 * ones(Nx, Ny);   % 4675 [m/s]
medium.density                 = 100 * ones(Nx, Ny);   % [kg/m^3]
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
t_end = 7e-6;
kgrid.makeTime(max(medium.sound_speed_compression(:),medium.sound_speed_shear(:)), cfl, t_end);
time = kgrid.t_array(kgrid.t_array<1*T_period);        %  5 cycle time array for sinusoidal input pulse
[t_sc, t_scale, t_prefix] = scaleSI(t_end);

%% Source Defination

% create initial pressure distribution using makeDis
%source.sxx = (1-2*(pi*source_f0*(kgrid.t_array - T_period)).^2 ).*exp(-(pi*source_f0.*(kgrid.t_array - T_period)).^2).*1/(pi*ds);
%source.sxx = (1-2*(pi^2*source_f0^2*(kgrid.t_array - T_period)).^2 ).*exp(-(pi^2*source_f0^2.*(kgrid.t_array - T_period)).^2);
%source.sxx = (1-2*(pi^2*source_f0^2*(kgrid.t_array - t_0)).^2 ).*exp(-(pi^2*source_f0^2.*(kgrid.t_array - t_0)).^2);
%source.sxx = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2).*sin(2*pi*source_f0*kgrid.t_array);
%source.syy = source.sxx;
%source.sxy = 0;
%disc1 = makeBall(Nx, Ny,Nz, Nx/2, Ny/2,4, 2); % define the sorce mask near the substrate surface

%% Pressure from temperature distribution

Pulse_duraton = 5e-9;
center = 10e-9;
laser_energy = 0.1e-3;
power = laser_energy/Pulse_duraton;
area = pi*(0.5e-3)^2;
intensity = power/area;


%  tem = temperature_time_gaussian_pulse(intensity_max,difusivity,conductivity,t0,tao,total_time)
% call tem.pressure(1,:)

% Diffusivity in m²/s and conductivity in W/m/K, 1cm-1 = 100 m-1, thermal
% conductivity from Scruby page 227

%tem = temperature_time_gaussian_pulse(1.6e10,0.88e-4,150,15e-9,10e-9,30e-9);     % For Si  
%tem =temperature_time_gaussian_pulse(1.6e10,1e-4,240,15e-9,10e-9,30e-9);       % For Al
%tem = temperature_time_gaussian_pulse(intensity,1.3e-5,50,10e-9,Pulse_duraton,70e-9);    % For 1% carbon Steel
tem = temp_init_pres(intensity,1.3e-5,50,0.5e-3,center,Pulse_duraton);
%tem = temp_init_pres(1.6e10,0.88e-4,150e3,1e-3,20e-9,10e-9);
press = squeeze(tem.pressure);
data_length = length(press(:,1));
source_mag = max(press(data_length/2,:));
%plot_time = kgrid.t_array(kgrid.t_array<= 30e-9);

%% Stress source DEfination

%
source.s_mask = zeros(Nx, Ny);
%source.s_mask= makeDisc(Nx,Ny, source_position_x, source_position_y, 5);
%source.s_mask = makeLine(Nx,Ny,[sub_grid_strt:sub_grid_strt+10,22],[sub_grid_strt:sub_grid_strt+10,32]);
%source.s_mask(sub_grid_strt:sub_grid_strt+1, 25:34) = 1;
source.s_mask(source_position_x, source_position_y:source_position_y+(data_length-1)) = 1;
%source.s_mask(:,:) = disc1 ;
%time = [0.33 , 0.57 ,2];

%
%source.sxx = source_mag * sin(2*pi*source_f0*time)
%max_pressure = max(tem.pressure(1,:));
%source.sxx = max_pressure*toneBurst(1/kgrid.dt, source_f0, 2);
%source.sxx = max_pressure*exp(-((kgrid.t_array - T_period)/(T_period/4)).^2);
%source.sxx = tem.pressure(1,:);
%source.sxx = 0.5*max_pressure;
%source.syy = source.sxx;
%source.sxy = 0;
%source_mag = max(source.sxx);
%}

%
for e = 1:data_length
    source.sxx(e,:) = press(e,2:end);
end
%}



%% Pressure Source deifination

%source.p0 =  makeDisc(Nx,Ny, 100, 30,3) ;
%source.p0(100,25:126) = source_mag * sin(2*pi*source_f0*time);
%source.p0 = source_mag * makeLine(Nx, Ny, [100, 30,],[100,40]);
%source.p0 = source_mag * makeDisc(Nx, Ny, 100, 30,3);

% same dimensions as the wedge
%gauss_space = exp(-((kgrid.x_vec-x0).^2 + (kgrid.y_vec-y0).^2)/ds);
%source.s_mask(:,:,1) = gauss_space;
%sensor.mask = makeCircle(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, 20);


%% Sensor Mask

%
%  Uncomment this line to apply the source for simulation

sensor.mask = zeros(Nx,Ny);

%{
source_sensor_distance = round(0.8*Ny);
%sensor.mask(100, source_sensor_distance) = 1;
sensor.mask(sub_grid_strt, source_sensor_distance) = 1;
sensor.mask(sub_grid_strt, round(0.4*Ny)) = 1;
source_sensor_dist = source_sensor_distance*dx;
%sensor.mask(100, 90) = 1;
%}

%
%sensors =  round(0.6*Ny):100:round(0.8*Ny);
%sensors =  2000:134:2600;
%sensors =  3010:200:3830;
%max_dist = Ny*dy;
%sensors =  round(31e-3/dx):round(2e-3/dx):round(33e-3/dx);
sensors = linspace(round(0.4*Ny),round(0.8*Ny),2);
source_sensor_dist = zeros(1,length(sensors));
i = 0; 
for k = sensors
%for k = slab_strt:1:slab_end
%    sensor.mask(slab_strt,k:k+10) = 1;               %2D
    sensor.mask(sub_grid_strt,k) = 1;               %2D
    %sensor.mask(k,round(0.8*Ny)) =1;
    i = i+1;
    source_sensor_dist(i) = (k - source_position_y)*dx;
end
%


%% define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
display_mask(sub_grid_strt, :) = 1;
%display_mask(crack_start_x:crack_end_x,crack_start_y:crack_end_y) = 1;         % Display the crack
%diameter = [5,10,15];
%source.p0 = max_pressure * makeDisc(Nx, Ny, sub_grid_strt+10, 30,5);
%source_mag = max_pressure;

%% Data Record
%sensor.record = {'u','u_split_field'};
sensor.record = {'p','u'};
% define input arguments
%
%input_args = {'PlotScale', [-1, 1] , 'PlotPML', false,...
%    'DisplayMask', display_mask,'MovieName','elastic_2d', 'DataCast', 'gpuArray-single','PMLAlpha',4};
input_args = {'PlotScale',[-source_mag, source_mag],'PlotPML', true, 'DisplayMask', display_mask + source.s_mask + sensor.mask, 'DataCast', 'gpuArray-single','PMLAlpha',4}; 
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
%{
input_args = {'PlotScale', [-source_mag, source_mag], 'PlotPML',false,...
              'DisplayMask', display_mask+source.s_mask+sensor.mask,...
              'RecordMovie', true, 'MovieName','surface_2d',...
              'DataCast','gpuArray-single'
              };
%}

%% run the simulation

sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
save('SRW_3points.mat');

%}
% =========================================================================
                               %% VISUALISATION
% =========================================================================
%
%[t_sc, t_scale, t_prefix] = scaleSI(t_end);
[p_sc, p_scale, p_prefix] = scaleSI(sensor_data.p(1,:));
[d_sc, d_scale, d_prefix] = scaleSI( source_sensor_dist(end));
    
%{
    %plot source, sensor, and position of the interface
    figure;
    imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, ...
        double(source.p0 | sensor.mask | display_mask), [-1 1]);
    colormap(getColorMap);
    ylabel('x-position [mm]');
    xlabel('y-position [mm]');
    axis image;
%}
%
% plot the re-ordered sensor data
%imagesc(sensor_data_reordered, [-1, 1]);
%subplot(2,1,1);

%{
% figure;
for k = 1:length(sensor_data.uz(:,1))
    %for k = 1:length(sensor_data.ux(:,1))
    %h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.uz(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
end
grid on;
legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Horizontal Velocity over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Velocity [m/s]');

%
figure;
%subplot(2,1,2);
for k = 1:length(sensor_data.ux(:,1))
    %for k = 1:length(sensor_data.ux(:,1))
    %h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
end
grid on;
legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'velocity Ux over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Velocity [m/s]');

figure;
%subplot(2,1,2);
for k = 1:length(sensor_data.uy(:,1))
    %for k = 1:length(sensor_data.ux(:,1))
    %h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.uy(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
end
grid on;
legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Velocity Uy over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Velocity [m/s]');

%{
    subplot(2,1,2);
    for k = 1:length(sensor_data.ux_split_s(:,1))
    %for k = 1:length(sensor_data.ux(:,1))
        %h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
        h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux_split_s(k,:));
        if k == 1, hold on, end
        legendinfo{k} = ['At sensor position = ' num2str(k)];
    end
%}
%}

%{
figure;
plot_no = length(sensor_data.p(:,1));
for k = 1:plot_no
    %for k = 1:length(sensor_data.ux(:,1))
    %h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.p(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
end
grid on;
legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Pressure over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [Pa]');
%}
%% For plotting a single time series data
%{
figure;
no_of_sensors_m = length(sensor_data.p(:,1));
for m = 1:no_of_sensors_m
    subplot(no_of_sensors_m,1,m);
    plot(kgrid.t_array * t_scale,sensor_data.p(m,:)*p_scale);
%    if m < no_of_sensors_m, set(gca,'XTick',[]), end
    b = source_sensor_dist*d_scale;
    grid on;
    legend( [' Source-Sensor distance = ' num2str(b)  d_prefix 'm'])
end
%set(gca,'xtickMode', 'auto')
%xticks(0:50:t_end*t_scale);
sgtitle('Pressure over Time') 
xlabel(['Time [' t_prefix 's]']);
ylabel(['Pressure [' p_prefix ' Pa]']);
grid on;
%}

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

% Plotting the stress input pulse
%{
figure;
plot(kgrid.t_array*t_scale,source.sxx);
xlabel(['Time [' t_prefix 's]']);
ylabel('Stress [Pa]');
title('Input pulse');
%plot(source.sxx);
grid on;
%}
%xticks(0:2:t_end*t_scale);

% Voxel plot
%voxelPlot(source.s_mask + sensor.mask);
%voxelPlot(p0 + sensor.mask);
%view([50, 50]);

%figure;
%voxelPlot(sensor.mask);
%title('Data record position');
%}
%surf();
%}

%{
figure;
[vec_sc, vec_scale, vec_prefix] = scaleSI(kgrid.x_vec);
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, ...
    double(source.s_mask | sensor.mask | display_mask), [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
%}

%% SOurce and sensor mask plot 
%{
figure;
[vec_sc, vec_scale, vec_prefix] = scaleSI(kgrid.x_vec);
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, ...
    double(source.s_mask | sensor.mask | display_mask), [-1 1]);
colormap(getColorMap);
ylabel(['x position [' vec_prefix 'm]']);
xlabel(['yposition [' vec_prefix 'm]']);
axis image;
%}