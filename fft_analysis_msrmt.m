clearvars
%msmnt_data_importing;
PA_Meas_Data_Demp;
signal_data = signal_data;
noof_rows = length(signal_data(:,1));

indx_strt_stl = 11000;
indx_end_stl = 40000;

indx_strt = indx_strt_stl;  % change the index based on the material
indx_end = indx_end_stl;
indx_move_stl =1800;

for r =1:noof_rows
    indx_strt = indx_strt + indx_move_stl; % chamge the indx_move for glass
    indx_end = indx_end + indx_move_stl;
    ms_data(r,:) = signal_data(r,indx_strt:indx_end);
    tim(r,:) = times(1,indx_strt:indx_end);
    a2(r,:) = fft(ms_data(r,:),2^15);
    dt = (tim(r,2)-tim(r,1));
    Fs = 1/dt;
    L= length(a2);
    f =Fs*(0:L/2)/L;
    freq_rng = 0.01*1e4;
    frequency = f(1:freq_rng);
    
    % P is not used yet
    p2(r,:) = abs(a2(r,:)/L);
    p1(r,:) = p2(r,1:L/2+1);
    p1(r,2:end-1) = 2*p1(r,2:end-1);
end

[t_sc, t_scale, t_prefix] = scaleSI(times(end));
[f_sc, f_scale, f_prefix] = scaleSI(frequency(end));
figure
%hold on
%
for k = 1:length(p1(:,1))
    subplot(noof_rows,1,k);
    plot(frequency*f_scale,p1(k,1:freq_rng));
  %  grid on;
end
%
xlabel(['Frequency [' f_prefix 'Hz]']);
ylabel(['Frequency Amplitude']);


%
figure
%hold on
for l = 1:length(ms_data(:,1))
    subplot(noof_rows,1,l);
    plot(tim(l,:)*t_scale,ms_data(l,:));
  %  grid on;
end
%}
