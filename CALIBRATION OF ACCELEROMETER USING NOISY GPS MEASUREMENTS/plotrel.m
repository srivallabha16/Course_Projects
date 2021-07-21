function plotrel(rel,tm,Edel,LP)
crel = randi(rel);
figure()
plot(tm, Edel(:,1,crel), 'k-', 'DisplayName', 'Position Error')
hold on
plot(tm, sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Position estimate')
title(sprintf('Error in Position estimate realization %d',crel))
legend
%saveas(gcf,sprintf('Error_in_position_estimate%d.png',crel))

figure()
plot(tm, Edel(:,2,crel), 'b-', 'DisplayName', 'Velocity Error')
hold on
plot(tm, sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Velocity estimate')
title(sprintf('Error in Velocity estimate realization %d',crel))
legend
%saveas(gcf,sprintf('Error_in_velocity_estimate%d.png',crel))

figure()
plot(tm, Edel(:,3,crel), 'm-', 'DisplayName', 'Bias Error')
hold on
plot(tm, sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in bias estimate')
title(sprintf('Error in bias estimate realization %d',crel))
legend
%saveas(gcf,sprintf('Error_in_bias_estimate%d.png',crel))
end