location =  'Lat69.945000_Lon-134.000000'; %'Lat65.004913_Lon40.054033'; %
foldername = 'SensitivityTestLocations'; 
%foldername = 'SuPerMapCircumArctic';

SuPerMap = load(strcat(foldername, filesep, 'SensitivityTest_DefaultValues', filesep, 'resultstruc_', location, '.mat'));
SuPerOCMap = load(strcat(foldername, filesep, 'OC_Run_oldversionstrata', filesep, 'Julia_', location, '.mat'));
Stratigraphy = load(strcat(foldername, filesep, 'strata_structs_oldversion', filesep, 'SuPerOCMap_', location, '.mat'));

T_Super = SuPerMap.resstruc.TS;
[depth_num, time_num] = size(T_Super);
depth = -SuPerMap.resstruc.pars(:,1);

time = SuPerOCMap.time;
Index_50k = length(time) - time_num + 1;

time_50k = time(Index_50k:end);
T_OC = SuPerOCMap.T(end-depth_num+1:end, Index_50k:end); %

% figure()
% for t=1:time_num
%     plot(T_Super(:,t), depth, T_OC(:,t), depth)
%     legend('SuPerMap', 'OC')
%     %pause(0.2)
%     drawnow
% end

disp('--- Lower Permafrost Boundary ---')
fprintf('SuPerMap Run: %.1f \n', depth(find(T_Super(:,end)<=0, 1, 'last')));
fprintf('OC Run: %.1f \n', depth(find(T_OC(:,end)<=0, 1, 'last')));

figure()
surf(T_Super, 'EdgeColor', 'none')
hold on
contour(T_Super, [0,0])
view(2)
set(gca, 'YDir', 'reverse')
title('SuPerMap')

figure()
surf(T_OC, 'EdgeColor', 'none')
hold on
contour(T_OC, [0,0])
view(2)
set(gca, 'YDir', 'reverse')
title('OC')

figure()
subplot(3,2,[1:4])
surf(time_50k, -depth, -abs(T_OC - T_Super), 'EdgeColor', 'none')
hold on
contour(time_50k, -depth, T_Super, [0,0], 'r')
contour(time_50k, -depth, T_OC, [0,0], 'b')
view(2)
set(gca, 'YDir', 'reverse')
title('Difference')
subplot(3,2,[5:6])
plot(time_50k, nanmean(T_OC - T_Super, 1))
axis tight


figure()
plot(SuPerMap.resstruc.timeF, SuPerMap.resstruc.TF, time, Stratigraphy.forcingData)
legend('SuPerMap', 'OC')