figure()
subplot(2,1,1)
plot(T_meas(1,:))
subplot(2,1,2)
plot(T_inv(1,:))


figure()
subplot(2,1,1)
imagesc(T_meas)
subplot(2,1,2)
imagesc(T_inv)