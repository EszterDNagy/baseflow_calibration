close all; clear variables; clc

% load sample dataset which contains one year of hourly precipitation 
% (p [mm]) and runoff (q [m^3/s]) data in a timetable
load('tt.mat')

% define a datetime array marking flood peaks which will be used for
% calibration
tpeak = datetime([2010 5 10 23 0 0; 2010 5 28 2 0 0; 2010 7 16 8 0 0; 2010 7 30 5 0 0; 2010 10 26 2 0 0]);

% define catchment area in km^2
A = 137;

% define the timestep of the time series in seconds
dt = 3600;

% define if the results should be plotted (1 - yes, 0 - no)
plot_mode = 1;

% run the calibration algorithm
[beta_0,Qb_0,Tc_0,Te_0,Tl_0,f] = beta_cal(tt.t,tt.q,tt.p,tpeak,A,dt,plot_mode);

% plot the ratio of Tc and Te to see if the values are between 1.5 and 3
figure
b = bar(Tc_0./Te_0,1);
b.FaceColor = [.5 .5 .5];
b.LineWidth = 1;
hold on
plot([.5 numel(Tc_0)+.5 NaN .5 numel(Tc_0)+.5],[1.5 1.5 NaN 3 3],'k--','linewidth',1.5)
plot([.5 numel(Tc_0)+.5],[2.25 2.25],'r-','linewidth',2)
ylim([0 3.2])
xlabel('Number of event')
ylabel('r = T_c/T_e [-]')
legend('r = T_c/T_e','accepted range','r_{opt} = 2.25','location','northoutside','orientation','horizontal')