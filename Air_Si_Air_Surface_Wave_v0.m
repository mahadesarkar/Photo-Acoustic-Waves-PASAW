% Simulation of elastic wave in 3D solid medium. 
clearvars;
% ========================================================================================================================================
%                                                                SIMULATION
% ========================================================================================================================================

% Properties of absorption layer
%}
pml_size = [10,10];
%lpha_y, alpha_z] = [20,20,0];
%}

%% Source frequency and time properties
source_mag = 5e3;
source_f0 = 10e6; % [Hz]
source_cycles = 2; % number of tone burst cycles
T_period = 1/(source_f0);
ppw       = 10;       % number of points per wavelength , use 1.5-3 for reyleigh wave
cfl   = 0.1;  
rayleigh_speed = 3000;

%% Grid parameters
dx = rayleigh_speed/ (ppw * source_f0);   % [m]
dy = dx;
dz = dx;
Nx = 64;           % number of grid points in the x (row) direction
Ny = 256;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% Medium Properties

%{
medium.sound_speed_compression = 343 * ones(Nx, Ny);   % 9138[m/s]
medium.sound_speed_shear       = 100 * ones(Nx, Ny);   % 4675 [m/s]
medium.density                 = 10 * ones(Nx, Ny);   % [kg/m^3]
%}

sub_grid_strt = round(0.3*Nx);      % Starting grid point of the Substrate medium
sub_grid_end = round(0.7*Nx);       % Starting grid point of the air medium
y_axis_sensor_pstn = round(0.1*Ny);


%
medium.sound_speed_compression = 930 * ones(Nx, Ny);   % 9138[m/s]
medium.sound_speed_shear       = 580 * ones(Nx, Ny);   % 4675 [m/s]
medium.density                 = 240 * ones(Nx, Ny);   % [kg/m^3]
%

medium.sound_speed_compression(sub_grid_strt:sub_grid_end,:) = 5960;   % [m/s]
medium.sound_speed_shear(sub_grid_strt:sub_grid_end,:)       = 3240;   % [m/s]
medium.density(sub_grid_strt:sub_grid_end,:)                 = 7900;   % [kg/m^3]

medium.sound_speed_compression(sub_grid_end:end,:) = 930;   % 9138[m/s]
medium.sound_speed_shear(sub_grid_end:end,:)       = 580;   % 4675 [m/s]
medium.density(sub_grid_end:end,:)                 = 240;   % [kg/m^3]
%}

% define the properties of the upper layer of the propagation medium
%{
medium.sound_speed_compression = 9300 * ones(Nx,Ny);   % [m/s]
medium.sound_speed_shear       = 5800 * ones(Nx,Ny);   % [m/s]
medium.density                 = 2398 * ones(Nx,Ny);   % [kg/m^3]
%}

%% define the absorption properties
medium.alpha_coeff_compression = 5 * ones(Nx,Ny) ; %0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 5 * ones(Nx,Ny);   %0.5; % [dB/(MHz^2 cm)]

medium.alpha_coeff_compression(sub_grid_strt:sub_grid_end,:) = 0.1 ; %0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear(sub_grid_strt:sub_grid_end,:)       = 0.1;   %0.5; % [dB/(MHz^2 cm)]
%medium.alpha_power = 2;      % Not erquired for Pstd_3d medium defination.
%medium.BonA =0.5;
%Nt = 1000;
%dt = 1e-5;

%% Time scaling 

%Nt = 1000;
%dt = 1e-5;
%kgrid.setTime(Nt, dt);
%t_end = 60*T_period;   % [s]
t_end = (Nx*dx)*10.5/rayleigh_speed;
kgrid.makeTime(max(medium.sound_speed_compression(:),medium.sound_speed_shear(:)), cfl, t_end);
time = kgrid.t_array(kgrid.t_array<1*T_period);        %  5 cycle time array for sinusoidal input pulse

% Scale the time and pressure values in SI unit
[t_sc, t_scale, t_prefix] = scaleSI(t_end);

%% Source Defination as stress source
%
%source.sxx = (1-2*(pi*source_f0*(kgrid.t_array - T_period)).^2 ).*exp(-(pi*source_f0.*(kgrid.t_array - T_period)).^2).*1/(pi*ds);
%source.sxx = (1-2*(pi^2*source_f0^2*(kgrid.t_array - T_period)).^2 ).*exp(-(pi^2*source_f0^2.*(kgrid.t_array - T_period)).^2);
%source.sxx = (1-2*(pi^2*source_f0^2*(kgrid.t_array - t_0)).^2 ).*exp(-(pi^2*source_f0^2.*(kgrid.t_array - t_0)).^2);
%source.sxx = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2).*sin(2*pi*source_f0*kgrid.t_array);
%source.syy = source.sxx;
%source.sxy = 0;
%disc1 = makeBall(Nx, Ny,Nz, Nx/2, Ny/2,4, 2); % define the sorce mask near the substrate surface

%slab = ones(Nx, Ny);
%slab (20:60,10:12)=1 ;
%source.s_mask = slab;

%{
source.s_mask = zeros(Nx, Ny);
source.s_mask= makeDisc(Nx,Ny, sub_grid_strt, 30,3);
source.sxx = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%source.s_mask= makeDisc(Nx,Ny, sub_grid_end, 30,3);
%source.s_mask = makeLine(Nx,Ny,[sub_grid_strt,25],[sub_grid_strt,45]);
%source.s_mask(:,:) = disc1 ;
%time = [0.33 , 0.57 ,2];
%}

%
source_mask = zeros(Nx,Ny);
%source_mask(sub_grid_strt :sub_grid_end,y_axis_sensor_pstn) = 1;
source_mask(sub_grid_strt,y_axis_sensor_pstn) = 1;
source_mask(sub_grid_end,y_axis_sensor_pstn) = 1;
source.s_mask = source_mask;
source.sxx(1,:) = source_mag*toneBurst(1/kgrid.dt, source_f0, source_cycles);
source.sxx(2,:) = -source_mag*toneBurst(1/kgrid.dt, source_f0, source_cycles);
%}

%source.sxx = source_mag;
%source.sxx = source_mag * sin(2*pi*source_f0*time);
%source.sxx(1,:) = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%source.sxx(2,:) = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%source.syy = source.sxx;
%source.sxy = 0;
%}

%{
source_mask=zeros(Nx,Ny);
source_mask(sub_grid_strt,26)=1;
source_mask(sub_grid_end,26)=1;
source.s_mask = source_mask;

% assign the source
source.sxx(1,:) = source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
source.sxx(2,:) = -source_mag*exp(-((kgrid.t_array - 2*T_period)/(T_period/2)).^2);
%}

%% Source as pressure source

%source.p0 =  makeDisc(Nx,Ny, 100, 30,3) ;
%source.p0(100,25:126) = source_mag * sin(2*pi*source_f0*time);
%source.p0 = source_mag * makeLine(Nx, Ny, [100, 30,],[100,40]);

%source.p0 = source_mag * makeDisc(Nx, Ny, sub_grid_strt, 30,4);

%% Sensor defination

% same dimensions as the wedge
%gauss_space = exp(-((kgrid.x_vec-x0).^2 + (kgrid.y_vec-y0).^2)/ds);
%source.s_mask(:,:,1) = gauss_space;
%sensor.mask = makeCircle(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, 20);
sensor.mask = zeros(Nx,Ny);
%source_sensor_distance1 = round(0.5*Ny);
%source_sensor_distance2 = round(0.6*Ny);
%source_sensor_distance3 = round(0.7*Ny);
%source_sensor_distance4 = round(0.8*Ny);
source_sensor_distance5 = round(0.9*Ny);

%sensor.mask(sub_grid_strt, source_sensor_distance1) = 1;
%sensor.mask(sub_grid_strt, source_sensor_distance2) = 1;
%sensor.mask(sub_grid_strt, source_sensor_distance3) = 1;
%sensor.mask(sub_grid_strt, source_sensor_distance4) = 1;
sensor.mask(sub_grid_strt, source_sensor_distance5) = 1;
%sensor.mask(source_sensor_distance, source_sensor_distance) = 1;

%% Define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
display_mask(sub_grid_strt, :) = 1;
display_mask(sub_grid_end, :) = 1;

%% Data record as pressure or velocity 

%sensor.record = {'u','u_split_field'};
sensor.record = {'p','u'};
%% Define input arguments
%{
input_args = {'PlotScale', [-1, 1] , 'PlotPML', false,'DisplayMask', display_mask,'MovieName','elastic_2d',...
    'MovieArgs', {'FrameRate', 10},'PlotFreq', 100,'MovieProfile', 'MPEG-4', 'DataCast', 'gpuArray-single','PMLAlpha',4};
%}

%
input_args = {  'PlotScale', [-source_mag, source_mag],...
                'PlotPML', false,...
                'DisplayMask', display_mask+source.s_mask+sensor.mask,...
                'DataCast', 'gpuArray-single',...
                'PMLSize', pml_size...
                'PMLAlpha',4    }; 
%}
%% run the simulation

sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
%}

% ==========================================================================================================================================
%                                                             VISUALISATION
% ==========================================================================================================================================

%% Plot the pressure values
%
figure;
no_of_sensors_k = length(sensor_data.p(:,1));
for k = 1:no_of_sensors_k
    subplot(no_of_sensors_k,1,k);
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.p(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
    grid on;
end

legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Pressure over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [Pa]');
%}

%% Plot the velocity in direction of x-axis
%
figure;
no_of_sensors_l = length(sensor_data.ux(:,1));
for l = 1:no_of_sensors_l
    subplot(no_of_sensors_l,1,k);
    h(k) = plot(kgrid.t_array * t_scale,sensor_data.ux(k,:));
    if k == 1, hold on, end
    legendinfo{k} = ['At sensor position = ' num2str(k)];
    grid on;
end

grid on;
legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Vertical Velocity Ux over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Velocity [m/s]');
%}

%% Plot the velocity in direction of y-axis
%
figure;
no_of_sensors_m = length(sensor_data.uy(:,1));
for m = 1:no_of_sensors_m
    subplot(no_of_sensors_m,1,m);
    h(m) = plot(kgrid.t_array * t_scale,sensor_data.uy(m,:));
    if m == 1, hold on, end
    legendinfo{m} = ['At sensor position = ' num2str(m)];
    grid on;
end

legend(h,legendinfo);
%xticks(0:2:t_end*t_scale);
title( 'Horizontal Velocity Uy over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel('Velocity [m/s]');

%}

%% Plot the Source and Sensor mask/position on the grid
%{
figure;
[vec_sc, vec_scale, vec_prefix] = scaleSI(kgrid.x_vec);
imagesc(kgrid.y_vec * vec_scale, kgrid.x_vec * vec_scale, ...
    double(source.s_mask     | sensor.mask), [-1 1]);
colormap(getColorMap);
ylabel(['x position [' vec_prefix 'm]']);
xlabel(['yposition [' vec_prefix 'm]']);
axis image;
%}
figure;
p=sensor_data.p;
subplot(2,2,1),imagesc(p),title('Preassure');
ux=sensor_data.ux;
subplot(2,2,2),imagesc(ux),title('Ux');
uy=sensor_data.uy;
subplot(2,2,3),imagesc(uy),title('Uy');