function display_Results(savename, plotdepth)
% to read and display temperatures from file

%in the moment, this interpolates the results to a "master grid"
%think about a better was to display results - maybe surf or patch?

if nargin <1
    savename = 'testlocation1_50k';
end
if nargin < 2
    plotdepth = -800;
end

load(strcat(savename, filesep, savename, '.mat'));


Time = RESULTS.time;
master_midpoints = RESULTS.depthInterp; %think about this. Needs to adapt somehow to the upperPos of the first module

TForcing = FORCING.TForcing;
saltConcForcing = FORCING.saltConcForcing;
timeForcing = FORCING.timeForcing;

T = RESULTS.T;
saltConc = RESULTS.saltConc;
liqWater = RESULTS.liqWater;
thermCond = RESULTS.thermCond;
c_eff = RESULTS.c_eff;

layerThick = master_midpoints(2:end) - master_midpoints(1:end-1);
saltConc_m = (saltConc(1:end-1,:) + saltConc(2:end,:) ) / 2;
liqWater_m = (liqWater(1:end-1,:) + liqWater(2:end,:) ) / 2;
totalSalt = sum(saltConc_m.*liqWater_m.*repmat(layerThick, 1, length(Time)), 1);




figure('Position', [1 1 1810 1227])
clf

subplot(7,1,2)
surf(Time, master_midpoints, T, 'EdgeColor', 'none')
view(2)
hold on
contour(Time, master_midpoints, T, [0,0], 'color', 'black')
colorbar
axis xy
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
title('Temperature over time')
ylabel('depth / m')

subplot(7,1,3)
surf(Time, master_midpoints, saltConc, 'EdgeColor', 'none')
view(2)
ylim([plotdepth,0])
set(gca, 'YDir', 'reverse')
colorbar
axis xy
xlim([timeForcing(1), timeForcing(end)])
title('Salt Concentration over time')
ylabel('depth / m')

%get axes position with colorbar
pos0 = get(gca, 'Position');

subplot(7,1,4)
plot(Time, totalSalt);
xlim([timeForcing(1), timeForcing(end)])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2) pos0(3) pos(4)]);
title('Total Salt over time')

subplot(7,1,5)
surf(Time, master_midpoints, liqWater, 'EdgeColor', 'none')
view(2)
ylim([plotdepth,0])
xlim([timeForcing(1), timeForcing(end)])
set(gca, 'YDir', 'reverse')
colorbar
axis xy
title('Liquid Water Content over time')
ylabel('depth / m')

subplot(7,1,6)
surf(Time, master_midpoints, thermCond, 'EdgeColor', 'none')
view(2)
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
title('Thermal conductivity over time')
colorbar
axis xy
ylabel('depth / m')

subplot(7,1,7)
surf(Time, master_midpoints, c_eff, 'EdgeColor', 'none')
view(2)
xlim([timeForcing(1), timeForcing(end)])
ylim([plotdepth,0])
axis xy
colorbar
caxis([0,4e7])
title('Heat capacity over time')
ylabel('depth / m')
xlabel('time / a')

ax1 = subplot(7,1,1);
%yyaxis(ax1,'left')
ax1 = plotyy(timeForcing, TForcing, timeForcing, saltConcForcing);
ylabel(ax1(1), 'temperature / ^\circ C')
%yyaxis(ax1,'right')
hold on
%plot(timeForcing, saltConcForcing)
ylabel(ax1(2), 'salt concentration / mol / m^3')
title('Forcing Data over time')
xlim([timeForcing(1), timeForcing(end)])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2) pos0(3) pos(4)]);


saveas(gcf, fullfile(savename, strcat(savename, '_julia_noSalt_results_overview.png')))


fprintf('Timesteps: %.0f \n', RUNINFO.timesteps);
%fprintf('Minimal Timestep: %3.2f days \n', RUNINFO.dt_min);
%fprintf('Maximum Timestep: %3.2f days \n', RUNINFO.dt_max);
fprintf('Runtime: %.2f hours \n', RUNINFO.endtime./(60*60));


% figure
% plot(TForcing)
% figure
% plot(D_surf')
% figure
% plot(density)

save(fullfile(savename, 'Results_Matrices.mat'), 'Time', 'master_midpoints', 'TForcing', 'saltConcForcing', ...
                                              'timeForcing', 'T', 'saltConc', 'thermCond', 'c_eff', ...
                                              'liqWater', 'totalSalt');