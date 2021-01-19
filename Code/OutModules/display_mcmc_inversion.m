close all
name = 'inversion_subsea_cutoffCosine_1_9';

results = load(sprintf('%s/%s.mat', name, name));

figure();clf

plot(results.measuredTemp,results.x_vec,'o', results.T_out, results.x_vec, '-')
%set(gca, 'YDir', 'reverse')
legend('data','dram estimate')

%cloud plots
figure();clf
plot(results.Chain(:,1),results.Chain(:,2),'.')
xlabel('p_1'); ylabel('p_2');title('MCMC chain');

try
    figure();clf
    plot(results.Chain(:,3),results.Chain(:,4),'.')
    xlabel('p_3'); ylabel('p_4');title('MCMC chain');
end

%chains
figure();clf
subplot(4,1,1)
plot(results.Chain(:,1),'.');ylabel('p_1');
title(sprintf('dram chain. Accepted %.1f%%',results.accepted*100))
%title(sprintf('\\tau = %.1f',tau(1)));
subplot(4,1,2)
plot(results.Chain(:,2),'.');ylabel('p_2');
%title(sprintf('\\tau = %.1f',tau(2)));
subplot(4,1,3)
plot(results.Chain(:,3),'.');ylabel('p_3');
%title(sprintf('\\tau = %.1f',tau(2)));
try
    subplot(4,1,4)
    plot(results.Chain(:,4),'.');ylabel('p_4');
    %title(sprintf('\\tau = %.1f',tau(2)));
end

%histograms
figure();clf
subplot(4,1,1)
histogram(results.Chain(:,1));
hold on
plot(results.p_exact(1), 0, 'ro', results.mean(1), 0 , 'yo'); ylabel('p_1')
subplot(4,1,2)
histogram(results.Chain(:,2));
hold on
plot(results.p_exact(2), 0, 'ro', results.mean(2), 0 , 'yo'); ylabel('p_2')
subplot(4,1,3)
histogram(results.Chain(:,3));
hold on
plot(results.p_exact(3), 0, 'ro', results.mean(3), 0 , 'yo'); ylabel('p_3')
try
    subplot(4,1,4)
    histogram(results.Chain(:,4));
    hold on
    plot(results.p_exact(4), 0, 'ro', results.mean(4), 0 , 'yo'); ylabel('p_4')
    legend('Ensemble', 'Exact Value', 'Reconstruction')
end


disp('--------------------------------------')
disp(strcat('--- ', 'dram', '---'))     
disp('--- Exact Main Parameter ---')
p_len = length(results.mean);
disp(results.p_exact)
disp('--- Reconstructed Main Parameter ---')
disp(results.mean)
disp('--- Final Reconstruction Error ---')
disp(norm(results.p_exact - results.mean)/norm(results.p_exact))
disp('--- Final Discrepancy ---')
disp(norm(results.measuredTemp - results.T_out)/norm(results.measuredTemp))