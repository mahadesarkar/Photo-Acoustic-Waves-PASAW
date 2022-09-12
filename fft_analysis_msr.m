clearvars
msmnt_data_importing;
noof_rows = length(signal_data(:,1));
indx_strt_Al = 1100;
indx_end_Al = 1900;
indx_strt_stl = 750;
indx_end_stl = 1200;
indx_strt_gls = 200;
indx_end_gls = 750;
indx_strt = indx_strt_stl;  % change the index based on the material
indx_end = indx_end_stl;
indx_move_gls = 80;
indx_move_Al =150;
indx_move_stl =100;

for r =1:noof_rows
    indx_strt = indx_strt + indx_move_stl; % chamge the indx_move for glass
    indx_end = indx_end + indx_move_stl;
    ms_data(r,:) = signal_data(r,indx_strt:indx_end);
    time(r,:) = times(1,indx_strt:indx_end);
    a2(r,:) = fft(ms_data(r,:),2^13);
    dt = (time(r,2)-time(r,1));
    Fs = 1/dt;
    L= length(a2);
    f =Fs*(0:L/2)/L;
    
    % P is not used yet
    p2(r,:) = abs(a2(r,:)/L);
    p1(r,:) = p2(r,1:L/2+1);
    p1(r,2:end-1) = 2*p1(r,2:end-1);
end

[t_sc, t_scale, t_prefix] = scaleSI(times(end));
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
    plot(time(l,:)*t_scale,ms_data(l,:));
  %  grid on;
end

