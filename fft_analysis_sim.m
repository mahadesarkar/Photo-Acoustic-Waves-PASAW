%clearvars;
%load('US_Al.mat');
signal_data = sensor_data.p;
times = kgrid.t_array;
noof_rows = length(signal_data(:,1));

indx_strt_Al = 10800;
indx_end_Al = 12500;

indx_strt_stl = 6700;
indx_end_stl = 9000;

indx_strt_gls = 8500;
indx_end_gls = 10050;

indx_strt_lsr = 3000;
indx_end_lsr = 7000;

indx_strt = indx_strt_gls;  % change the index based on the material
indx_end = indx_end_gls;

indx_move_gls = 800;
indx_move_Al = 950;
indx_move_stl = 950;
indx_move_lsr = 1550;

clear ms_data;
clear tim;
for r =1:noof_rows
    indx_strt = indx_strt + indx_move_gls; % chamge the indx_move for glass
    indx_end = indx_end + indx_move_gls;
    ms_data(r,:) = gather(signal_data(r,indx_strt:indx_end));
    tim(r,:) = times(1,indx_strt:indx_end);
    a2(r,:) = fft(ms_data(r,:),2^14);
    dt = (tim(r,2)-tim(r,1));
    Fs = 1/dt;
    L= length(a2);
    f =Fs*(0:L/2)/L;
    
    % P is not used yet
    p2(r,:) = abs(a2(r,:)/L);
    p1(r,:) = p2(r,1:L/2+1);
    p1(r,2:end-1) = 2*p1(r,2:end-1);
end

[f_sc, f_scale, f_prefix] = scaleSI(f(end));
figure
%hold on
for k = 1:length(p1(:,1))
    subplot(noof_rows,1,k);
    plot(f*f_scale,p1(k,:));
  %  grid on;
end

xlabel(['Frequency [' f_prefix 'Hz]']);
ylabel(['Frequency Amplitude']);


figure
%hold on
for l = 1:length(ms_data(:,1))
    subplot(noof_rows,1,l);
    plot(tim(l,:)*t_scale,ms_data(l,:));
  %  grid on;
end




[p_sc, p_scale, p_prefix] = scaleSI(sensor_data.p(1,:));
[d_sc, d_scale, d_prefix] = scaleSI( source_sensor_dist(end));

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

