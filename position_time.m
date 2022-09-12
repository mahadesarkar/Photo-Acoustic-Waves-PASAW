clearvars
%positions = [0 10 20 30 40 50 60 70 80 90 100 110]*1E-3; % m
%peak_times = [32.100 34.890 39.030 41.720 45.040 49.600 53.220 56.140 59.300 62.240 65.370 68.610]*1E-6;
%positions = ([0 10 20 30 40 50 60 70 80 90 ])*1E-3; % m
positions = ([0 2 4 6 8 10 12 14 16 18 ])*1E-3; % m
%positions = [0 5 10 15 20 25 30 35 40 45 ]*1E-3; % m
%peak_times_Al_measured = [35.3 38.1 41.6 45.7 48.1 51.74 56.24 58.76 61.73 64.85]*1E-6;
%peak_times_steel_measured = [26.7 28.5 30.4 32 33.64 35.32 37.1 38.84 40.54 42.24]*1E-6;
%peak_times_glass_measured = [19.2 20.52 22 23.68 25.26 26.8 28.48 30.08 31.64 33.24]*1E-6;
%peak_times = [peak_times_Al;peak_times_steel;peak_times_glass];
%peak_times_Al_simulated = [31.56 35.12 38.6 42.06 45.5 49.04 52 55.1 59.1 62.4]*1E-6;
%peak_times_steel_simulated = [31.22 34.93 38.35 41.78 45.17 48.46 51.86 55.28 58.64 62.15]*1E-6;
%peak_times_glass_simulated = [31.41 34.76 38.09 41.45 44.6 47.58 51.2 54.4 57.07 60.1]*1E-6;
peak_times_Al_simulated_lsr = [29.62 33.06 36.6 39.65 43.11 46.76 49.83 53.73 57.17 60.38]*1E-6;
simulated_gls = [1, 1.267, 1.20, 1.32, 1.45, 1.6, 1.72, 1.77, 1.76, 1.93]*1e-6;
simulated_Al = [1.19, 1.03, 1.30,1.45, 1.47, 1.63, 1.72, 1.83, 1.92,1.96]*1e-6;
simulated_stl = [1.76, 1.76, 1.66, 1.61,1.54,1.47,1.51, 1.42,1.42, 1.39]*1e-6;
simulated_stl_lsr = [8.32, 8.96, 9.69,10.31,11.06,11.68, 12.40,13.06,13.77,14.40]*1e-6;
peak_frq_stl_lsr_sim = [2, 1.89,1.85,1.83, 1.8, 1.75, 1.71, 1.7, 1.69, 1.69 ]*1e6;
peak_frq_stl_lsr_sim = [2, 1.89,1.85,1.83, 1.8, 1.75, 1.71, 1.7, 1.69, 1.69 ]*1e6;

%
%[t_sc, t_scale, t_prefix] = scaleSI(peak_times_Al_simulated_lsr(1));
[t_sc, t_scale, t_prefix] = scaleSI(simulated_stl_lsr(1));
[vec_sc, vec_scale, vec_prefix] = scaleSI(positions);

%[xData, yData] = prepareCurveData(peak_times_Al_simulated_lsr*t_scale, positions*vec_scale );
[xData, yData] = prepareCurveData(simulated_stl_lsr*t_scale, positions*vec_scale );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );


f2 = figure;
h = plot( fitresult, xData, yData );
grid on;
legend( h, 'positions vs. peak_times', 'lin. regr.', 'Location', 'best', 'Interpreter', 'none' );
annotation(f2,'textbox',...
    [0.144214285714286 0.743333333333333 0.295785714285714 0.163333333333337],...
    'String',{'f(x) = p1*x + p2','Coefficients (95%cb):','p1 = 2993  (2839, 3148)','p2 = -0.1048(-0.1127, -0.09693)','Goodness of fit: ','SSE: 3.292e-05',' R-square: 0.996','Adjuste R-square: 0.9955','RMSE: 0.002029'},...
    'FitBoxToText','off',...
    'FontSize',8);
    
xlabel( 'peak_times', 'Interpreter', 'none' );
ylabel( 'positions', 'Interpreter', 'none' );

xlabel( ['Time [' t_prefix 's]']);
ylabel( ['Position [' vec_prefix 'm]']);
velocity = fitresult.p1*1e3 ;
%title(['Rayleigh Wave Velocity ', num2str(velocity) , 'm/s']);
title(['Wave Velocity ', num2str(round(velocity)) , 'm/s']);

%set(gca, 'XAxisLocation', 'origin');
grid on


%{
figure;
for i=1:3

    [t_sc, t_scale, t_prefix] = scaleSI(peak_times(i,1));
    [vec_sc, vec_scale, vec_prefix] = scaleSI(positions);
    
    [xData, yData] = prepareCurveData( peak_times(i,:)*t_scale, positions*vec_scale );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    
    
    h = plot( fitresult, xData, yData );
    grid on;
    legend( h, 'positions vs. peak_times', 'lin. regr.', 'Location', 'best', 'Interpreter', 'none' );
   % annotation(figure,'textbox',...
   %     [0.144214285714286 0.743333333333333 0.295785714285714 0.163333333333337],...
   %     'String',{'f(x) = p1*x + p2','Coefficients (95%cb):','p1 = 2993  (2839, 3148)','p2 = -0.1048(-0.1127, -0.09693)','Goodness of fit: ','SSE: 3.292e-05',' R-square: 0.996','Adjuste R-square: 0.9955','RMSE: 0.002029'},...
   %     'FitBoxToText','off',...
   %     'FontSize',8);
        
        
    xlabel( 'peak_times', 'Interpreter', 'none' );
    ylabel( 'positions', 'Interpreter', 'none' );
    
    xlabel( ['Time [' t_prefix 's]']);
    ylabel( ['Position [' vec_prefix 'm]']);
    velocity = fitresult.p1*1e3 ;
    %title(['Rayleigh Wave Velocity ', num2str(velocity) , 'm/s']);
    title(['Surface Wave Velocity ', num2str(round(velocity)) , 'm/s']);
    grid on
    hold on
end
%}