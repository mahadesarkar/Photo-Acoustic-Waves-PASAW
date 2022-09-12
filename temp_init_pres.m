% For the simulation of equation 3.7 from book Ready and my printout from
% professor's book C.B. Scruby eqn 5.27, Rectangular wave pulse laser temperature generation.
% Call the function tem = temperature_time_v5(1.9e10,0.88e-4,148,0.3e-3,20e-9)
% temperature_time_gauss_3d(1.6e10,0.88e-4,150e3,1e-3,20e-9,10e-9)
% tem = temp_init_pres(intensity,1.3e-5,50,0.5e-3,center,Pulse_duraton);
% With the material properties of steel

function tem = temp_init_pres(intensity_max,difusivity,conductivity,beam_radius,delay,tao)

%
thermal_expnsion = 2.5e-6;  % for Si3N4 3.5e-6/k, SiO2 0.5e-6/K, Si 2.5, Ge 5.8
                            % https://www.tf.uni-kiel.de/matwis/amat/semi_en/kap_3/backbone/r3_1_3.html
compressibility = 4.8e-12;  % k =1/(rho*velocity²)


% Grid setting
% Grid size should be same as Acoustic grid size
rayleigh_speed = 3000;
ppw = 7;
source_f0 = 2e6;
Nx = 64;            % number of grid points in the x (row) direction
%Ny = Nx;           % number of grid points in the y (column) direction
Ny = 4;            % Number of points on the line source and determine from the laser spot diameter
Nz = Nx/4;          % Less number of grid points, need onlz the first or first few points
dx = rayleigh_speed/ (ppw * source_f0);
%dx = 0.1e-3;        % grid point spacing in the x direction [m] 0.1mm
dy = dx;            % grid point spacing in the y direction [m] 0.1mm
%dz = 0.5e-6;
dz = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);


%medium.sound_speed = 6000*ones(Nx,Ny);
cfl   = 0.1;    % Courant-Friedrichs-Lewy number
t_end = 70e-9;

kgrid.makeTime(rayleigh_speed, cfl, t_end);


%% assign only the necessary grid values  

tem.x = kgrid.x_vec(Nx/2+1); % For the line source choosing onle the center line of the laser spot
tem.y = kgrid.y_vec(:);      % all values in y-axis as line source is along y-axis
tem.z = kgrid.z_vec(9);      % only zero for the surface temperature
tem.time_axis = kgrid.t_array(2:end)';  % xclude the first value as it is zero which will cause nan in calculation
tem.intensity = intensity_max;
tem.beam_radius = beam_radius;
tem.tao = tao;
tem.t0 = delay;

tem.temperature = zeros(length(tem.x),length(tem.y),length(tem.z),length(tem.time_axis)); %  3d temperature variable
tem.pressure = zeros(length(tem.x),length(tem.y),length(tem.z),length(tem.time_axis));
tem.gaus = zeros(length(tem.x),length(tem.y),length(tem.z));
pressure_prefactr = thermal_expnsion/compressibility;

tem.prefactor = (tem.intensity*sqrt(difusivity/pi)*tem.beam_radius^2)/conductivity;

for k = 1:length(tem.z)
    for i = 1:length(tem.x)
        for n = 1:length(tem.y)
 
        fun = @(t_bar,t) exp(-2.77.*((t -tem.t0 - t_bar )/tem.tao).^2).*exp(((-(tem.z(k)^2))./(4*difusivity.*t_bar)) - ((tem.x(i)^2 + tem.y(n)^2)./(4*difusivity.*t_bar + tem.beam_radius^2))).*(1./(sqrt(t_bar).*(4*difusivity.*t_bar+tem.beam_radius^2)));
        %fun = @(t_bar,t) exp(((-(tem.z(k)^2))./(4*difusivity.*t_bar)) - ((tem.x(i)^2 + tem.y(n)^2)./(4*difusivity.*t_bar + tem.beam_radius^2))).*(1./(sqrt(t_bar).*(4*difusivity.*t_bar+tem.beam_radius^2)));
            for time1 = 1:length(tem.time_axis)
                tem.gaus(i,n,k) = exp(-(tem.x(i)^2 + tem.y(n)^2)./(tem.beam_radius^2));
                t = tem.time_axis(time1);
                tem.temperature(i,n,k,time1) = tem.prefactor*integral(@(t_bar) fun(t_bar,t), 0, t);
                tem.pressure(i,n,k,time1) = pressure_prefactr*tem.temperature(i,n,k,time1);
            end      
        end
    end
end

end

%pressur = squeeze(tem.pressure);
%{
tem.temperature(:,:,:,1) = 0;
[t_sc, t_scale, t_prefix] = scaleSI(tem.time_axis);
[tem_sc, tem_scale, tem_prefix]= scaleSI(squeeze(tem.temperature(1,1,1,:)));
h = zeros(size(tem.z));
for w = 1:length(tem.z)
    %h(w) = plot(tem.time_axis*t_scale,squeeze(tem.temperature(36,36,w,:)));
    h(w) = plot(tem.time_axis*t_scale,squeeze(tem.temperature(1,1,w,:))*tem_scale);
    %legend(legendinfo)
    if w == 1, hold on, end
    [x_sc, scale, prefix] = scaleSI(tem.z(w));
    legendinfo{w} = ['At depth Z = ' num2str(x_sc), 'm'];
    grid on;
    %xticks(0:10:max(tem.time_axis));
    %yticks(0:25:325);
    
end


legend(h,legendinfo);
title( 'Temperature response over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel(['Temperature rise / [' tem_prefix 'K]']);
figure;
%surf(squeeze(tem.gaus(:,:,1)));
%title("Gaussian Input Pulse")
%}

%{
figure;
tem.pressure(:,:,:,1) = 0;
[pre_sc, pre_scale, pre_prefix]= scaleSI(squeeze(tem.pressure(1,1,1,:)));
h = zeros(size(tem.z));
for w = 1:length(tem.z)
    %h(w) = plot(tem.time_axis*t_scale,squeeze(tem.temperature(36,36,w,:)));
    h(w) = plot(tem.time_axis*t_scale,squeeze(tem.pressure(1,1,w,:))*pre_scale);
    %legend(legendinfo)
    if w == 1, hold on, end
    [x_sc, scale, prefix] = scaleSI(tem.z(w));
    legendinfo{w} = ['At depth Z = ' num2str(x_sc), 'm'];
    grid on;
    %xticks(0:10:max(tem.time_axis));
    %yticks(0:25:325);
    
end


legend(h,legendinfo);
title( 'Pressure response over Time' );
xlabel(['Time [' t_prefix 's]']);
ylabel(['Pressure rise / [' pre_prefix 'Pa]']);
figure;
%surf(squeeze(tem.gaus(:,:,1)));
%title("Pressure response over time")

%}

%{

[t_sc, t_scale, t_prefix] = scaleSI(t_end);
figure;
plot(kgrid.t_array * t_scale,tem.sensor_data);
%plot(kgrid.t_array * t_scale,tem.sensor_data);

grid on;
%}