clearvars
close all

om = 0.1;
ts = 1/200;% sampling frequency of 200Hz
t = 0:ts:30;
tm = 0:1/5:30;

wm = 0; ws = 0.02; %accelerometer noise mean and sigma
bm = 0; bs = 0.1; %bias mean and sigma
pm = 0; ps = 10; %position meand and sigma
vm = 100; vs = 1; %velocity mean and sigma
pnm = 0; pns = 1; %x noise mean and sigma
vnm = 0; vns = 0.04; %v noise mean and sigma

rel = 1e4; %no of realizations
wn = zeros(rel, length(t));

b = normrnd(bm, bs, [rel,1]);
pini = normrnd(pm, ps,[rel,1]);
vini = normrnd(vm, vs,[rel,1]);

pn = zeros(rel, length(tm));
vn = zeros(rel, length(tm));

for tind = 1:length(t)
    wn(:, tind) = normrnd(wm, ws, [rel,1]);
end
histogram(wn(:,1))
for tind = 1:length(tm)
    pn(:,tind) = normrnd(pnm, pns, [rel,1]);
    vn(:,tind) = normrnd(vnm, vns, [rel,1]);
end

zk = zeros(2,1);

LXk = zeros(length(tm),3, rel);
LX = zeros(length(tm),3, rel);
Lino = zeros(length(tm),2, rel);
Edel = zeros(length(tm),3,rel);
LP = zeros(3,3,length(tm));

%Kalman Filter: Error variance propagation.
phi = [1 ts -0.5*ts*ts; 0 1 -ts;0 0 1];
gam = [-0.5*ts^2; -ts; 0];
H = [1 0 0; 0 1 0];
Xprop = [0;0;0];
Mprop = [ps^2 0 0;0 vs^2 0;0 0 bs^2];
V = [pns^2 0;0 vns^2];

for tind = 1:length(tm)
    K = Mprop * H' * inv(H * Mprop * H' + V);
    P = (eye(3) - K * H) * Mprop * (eye(3) - K * H)' + K * V * K';
    LP(:,:,tind) = P;
    
    Mprop = phi * P * phi' + gam * ws^2 * gam';
    for ind = 1:39
        Mprop = phi * Mprop * phi' + gam * ws^2 * gam';
    end
end
figure()
plot(squeeze(LP(1,1,:)))
title('Kalman P (1,1) element')

a = 10 * sin(om * t);
for ri = 1:rel
    p = zeros(length(t),1);
    v = zeros(length(t),1);
    ac = zeros(length(t),1);
    vc = zeros(length(t),1);
    pc = zeros(length(t),1);
    zp = zeros(length(tm),1);
    zv = zeros(length(tm),1);
    
    p(1) = pini(ri);
    v(1) = vini(ri);
    
    vc(1) = vm;
    pc(1) = pm;
    ac(1) = a(1) + b(ri) + wn(ri, 1);
    
  
    for tind = 2:length(t)
        v(tind) = v(1) + 10/om - 10/om * cos(om*t(tind));
        p(tind) = p(1) + (v(1) + 10/om) * t(tind) - 10/(om^2) * sin(om*t(tind));
        
        ac(tind) = a(tind) + b(ri) + wn(ri, tind);
        vc(tind) = vc(tind-1) + ac(tind-1) * ts;
        pc(tind) = pc(tind-1) + vc(tind-1)* ts + 0.5 * ac(tind-1) * ts^2;
    end
    
    zv = v(1:40:length(t),1) - vc(1:40:length(t),1) + (vn(ri,:))';
    zp = p(1:40:length(t),1) - pc(1:40:length(t),1) + (pn(ri,:))';
          
    %Kalman filter:finding the estimate
    Xk = [0;0;0];
    %Lpk(ri, 1) = Xk(1); Lvk(ri, 1) = Xk(2); Lbk(ri, 1) = Xk(3);
    LXk(1,:, ri) = Xk;
    
    for tind = 1:length(tm)
        
        zk = [zp(tind);zv(tind)];%IMPORTANT: Measurement step.
        
        Xa = [p(1+(tind-1)*40)-pc(1+(tind-1)*40); v(1+(tind-1)*40)-vc(1+(tind-1)*40); b(ri)];
        Lino(tind,:,ri) = zk - H * Xk;
        LX(tind,:,ri) = Xa;
        
        
        %Eino(tind,:,:) = squeeze(Eino(tind,:,:)) + (zk - H * Xk)*(zk - H * Xk)';
        
        Xk = Xk + squeeze(LP(:,:,tind)) * H' * inv(V) * (zk - H * Xk);
        LXk(tind,:,ri) = Xk;
        Edel(tind,:,ri) = Xa - Xk;
        Xk = phi^40 * Xk; %IMPORTANT: Propagation step.            

    end
    
end

%Avg error in estimate, across simulations
Er_ave = zeros(length(tm),3);
for tind = 1:length(tm)
    for ri = 1:rel
        Er_ave(tind,:) = (Er_ave(tind,:))' + (Edel(tind,:,ri))';
    end
end
Er_ave = Er_ave./rel;

P_ave = zeros(3,3,length(tm));
OErEs = zeros(3,3,length(tm));
for tind = 1:length(tm)
    for ri = 1:rel
        dd = (Edel(tind,:,ri))' - (Er_ave(tind,:))';
        P_ave(:,:,tind) = squeeze(P_ave(:,:,tind)) + dd * dd';
        OErEs(:,:,tind) = squeeze(OErEs(:,:,tind)) + dd * (LXk(tind,:,ri));
    end
end
P_ave = P_ave./rel;
OErEs = OErEs./rel;

Ores = zeros(2,2);
tindm = 30; tindi = 120;
for ri = 1:length(rel)
    Ores = Ores + (Lino(tindi,:,ri))' * Lino(tindm,:,ri);
end
Ores = Ores./rel;
disp('Ores')
disp(Ores)


%{
crel = randi(rel);
figure()
plot(tm, Edel(:,1,crel), 'k-', 'DisplayName', 'Position Error')
hold on
plot(tm, sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Position estimate')
title('Error in Position estimate realization 1')
legend
saveas(gcf,'Error_in_position_estimate1.png')

figure()
plot(tm, Edel(:,2,crel), 'b-', 'DisplayName', 'Velocity Error')
hold on
plot(tm, sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Velocity estimate')
title('Error in Velocity estimate realization 1')
legend
saveas(gcf,'Error_in_velocity_estimate1.png')

figure()
plot(tm, Edel(:,3,crel), 'm-', 'DisplayName', 'Bias Error')
hold on
plot(tm, sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in bias estimate')
title('Error in bias estimate realization 1')
legend
saveas(gcf,'Error_in_bias_estimate1.png')
%}
plotrel(rel,tm,Edel,LP);

figure()
subplot(3,1,1)
plot(tm, Er_ave(:,1),'r-','DisplayName','Avg error in position')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
title('Ensemble average of errors')
legend
subplot(3,1,2)
plot(tm, Er_ave(:,2),'b-','DisplayName','Avg error in velocity')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
subplot(3,1,3)
plot(tm, Er_ave(:,3),'m-','DisplayName','Avg error in bias')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
saveas(gcf, 'Simulated_Average_error.png')


figure()
subplot(3,1,1)
plot(tm,squeeze(P_ave(1,1,:)), 'b-', 'DisplayName','(1,1) Simul')
hold on
plot(tm, squeeze(LP(1,1,:)), 'r.','DisplayName','(1,1) Kalman')
title('Simulated vs Kalman Variance (diagonal terms)')
legend
subplot(3,1,2)
plot(tm,squeeze(P_ave(2,2,:)), 'b-',  'DisplayName','(2,2) Simul')
hold on
plot(tm, squeeze(LP(2,2,:)), 'r.','DisplayName','(2,2) Kalman')
legend
subplot(3,1,3)
plot(tm,squeeze(P_ave(3,3,:)), 'b-', 'DisplayName','(3,3) Simul')
hold on
plot(tm, squeeze(LP(3,3,:)), 'r.','DisplayName','(3,3) Kalman')
legend
saveas(gcf,'Simulated_Variance_diagonal.png')

figure()
subplot(3,1,1)
plot(tm,squeeze(P_ave(1,2,:)), 'b-', 'DisplayName','(1,2) Simul')
hold on
plot(tm, squeeze(LP(1,2,:)), 'r.','DisplayName','(1,2) Kalman')
title('Simulated vs Kalman Variance (off-diagonal terms)')
legend
subplot(3,1,2)
plot(tm,squeeze(P_ave(1,3,:)), 'b-', 'DisplayName','(1,3) Simul')
hold on
plot(tm, squeeze(LP(1,3,:)), 'r.','DisplayName','(1,3) Kal')
legend
subplot(3,1,3)
plot(tm,squeeze(P_ave(2,3,:)), 'b-', 'DisplayName','(2,3) Simul')
hold on
plot(tm, squeeze(LP(2,3,:)), 'r.','DisplayName','(2,3) Kal')
legend
saveas(gcf,'Simulated_Variance_off-diagonal.png')

figure()
subplot(4,1,1)
plot(tm,squeeze(OErEs(1,1,:)), 'r.','DisplayName','(1,1)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
title('Orthogonality of error & estimates')
legend
subplot(4,1,2)
plot(tm,squeeze(OErEs(1,2,:)), 'r.','DisplayName','(1,2)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
subplot(4,1,3)
plot(tm,squeeze(OErEs(1,3,:)), 'r.','DisplayName','(1,3)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
subplot(4,1,4)
plot(tm,squeeze(OErEs(2,3,:)), 'r.','DisplayName','(2,3)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
saveas(gcf, 'Ortho_Errors_Estimates.png')

figure()
subplot(4,1,1)
plot(tm,squeeze(OErEs(2,2,:)), 'r.','DisplayName','(2,2)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
title('Orthogonality of error & Estimates contd.')
legend
subplot(4,1,2)
plot(tm,squeeze(OErEs(2,1,:)), 'r.','DisplayName','(2,1)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
subplot(4,1,3)
plot(tm,squeeze(OErEs(3,1,:)), 'r.','DisplayName','(3,1)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
subplot(4,1,4)
plot(tm,squeeze(OErEs(3,3,:)), 'r.','DisplayName','(3,3)')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
legend
saveas(gcf, 'Ortho_Errors_Estimates_contd.png')

figure()
plot(tm, sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '+ 1 sigma')
hold on
plot(tm, 0*ones(length(tm),1),'k-','DisplayName','zero')
plot(tm, -sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '- 1 sigma')

plot(tm, Er_ave(:,2),'g.','DisplayName', 'Simulated Avg error in velocity estimate')
legend
title('Velocity and its errors')
saveas(gcf, 'Velocity_errors.png')

figure()
crel = randi(rel);
plot(tm, squeeze(LX(:,2,crel)),'b.', 'DisplayName','Simulated Velocity')
hold on
plot(tm, squeeze(LXk(:,2,crel)),'b-', 'DisplayName','Estimate of Velocity')
legend
title('Estimate and simulated velocity for a realisation')
saveas(gcf,'Est_and_simul_velo.png')

%{
crel = randi(rel);
figure()
plot(tm, Edel(:,1,crel), 'k-', 'DisplayName', 'Position Error')
hold on
plot(tm, sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Position estimate')
title('Error in Position estimate realization 2')
legend
saveas(gcf,'Error_in_position_estimate2.png')

figure()
plot(tm, Edel(:,2,crel), 'b-', 'DisplayName', 'Velocity Error')
hold on
plot(tm, sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Velocity estimate')
title('Error in Velocity estimate realization 2')
legend
saveas(gcf,'Error_in_velocity_estimate2.png')

figure()
plot(tm, Edel(:,3,crel), 'm-', 'DisplayName', 'Bias Error')
hold on
plot(tm, sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in bias estimate')
title('Error in bias estimate realization 2')
legend
saveas(gcf,'Error_in_bias_estimate2.png')

crel = randi(rel);
figure()
plot(tm, Edel(:,1,crel), 'k-', 'DisplayName', 'Position Error')
hold on
plot(tm, sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(1,1,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Position estimate')
title('Error in Position estimate realization 3')
legend
saveas(gcf,'Error_in_position_estimate3.png')

figure()
plot(tm, Edel(:,2,crel), 'b-', 'DisplayName', 'Velocity Error')
hold on
plot(tm, sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(2,2,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in Velocity estimate')
title('Error in Velocity estimate realization 3')
legend
saveas(gcf,'Error_in_velocity_estimate3.png')

figure()
plot(tm, Edel(:,3,crel), 'm-', 'DisplayName', 'Bias Error')
hold on
plot(tm, sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '+ 1 sigma')
plot(tm, -sqrt(squeeze(LP(3,3,:))), 'r-', 'DisplayName', '- 1 sigma')
xlabel('Measurement time')
ylabel('Error in bias estimate')
title('Error in bias estimate realization 3')
legend
saveas(gcf,'Error_in_bias_estimate3.png')
%}


