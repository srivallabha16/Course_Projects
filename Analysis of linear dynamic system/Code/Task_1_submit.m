%****************  TASK 1   *********************
clear
close all

load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs

load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);

%****************  TASK 1, Question 1 and 2  *********************

ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;

m = 2; %number of output channels
q = 2; %number of input channels
% n s.t size(H) = mn X qn
n = 100;
%n = (401 - 41)/2 ;%n is half the number of samples beginning from t = 1

ind_off = 41;
H = zeros(m*n, q*n);
Htil = zeros(m*n, q*n);

%To create H and H tilda matrices
r = 1;
c = 1;
for k = 1 : n
    ind = ind_off + k + [0:n-1];     
    
    H(r,[c:q:c+q*(n-1)]) = y11(ind);
    H(r+1,[c:q:c+q*(n-1)]) = y21(ind);
    H(r,[c+1:q:c+1+q*(n-1)]) = y12(ind);
    H(r+1,[c+1:q:c+1+q*(n-1)]) = y22(ind);
    
    Htil(r,[c:q:c+q*(n-1)]) = y11(ind+1);
    Htil(r+1,[c:q:c+q*(n-1)]) = y21(ind+1);
    Htil(r,[c+1:q:c+1+q*(n-1)]) = y12(ind+1);
    Htil(r+1,[c+1:q:c+1+q*(n-1)]) = y22(ind+1);
    
    r = r + m;
end

s_vals = svd(H);
figure(1);
plot(s_vals,'b.','MarkerSize',10)
xlabel('singular value index','FontSize',12);
ylabel('H_{100} singular values','FontSize',12);
grid on

%when ns = 6
ns = 6;
[U, S, V] = svd(H);
U1 = U([1:m*n],[1:ns]);
S1 = S([1:ns],[1:ns]);
V1 = V([1:n*q],[1:ns]);

On6 = U1 * sqrt(S1);
Cn6 = sqrt(S1) * V1';

On6inv = inv(sqrt(S1)) * U1';
Cn6inv = V1 * inv(sqrt(S1));

C6 = On6([1:m],:);
B6 = Cn6(:,[1:q]);
A6 = On6inv * Htil * Cn6inv;

mev6 = max(abs(eig(A6)));
fprintf('When ns = 6, dominant eval is: %f \n',mev6);

%Plotting impulse response
h = zeros(m,q,n);
x = B6;
for k = 1:n
    h(:,:,k) = C6 * x;
    x = A6 * x;
end

tsim = [1 : n] * ts;
figure()
subplot(211)
plot(tsim,squeeze(h(1,1,:)),'b*', t, y11, 'g*')
xlabel('Time')
ylabel('y_{11}')
legend('Model6', 'Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,1,:)),'b*', t, y21, 'g*')
xlabel('Time')
ylabel('y_{21}')
legend('Model6','Data')
grid on
axis([-0.2 2 -0.1 0.1])

figure()
subplot(211)
plot(tsim,squeeze(h(1,2,:)),'b*', t, y12, 'g*')
xlabel('Time')
ylabel('y_{12}')
legend('Model6','Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,2,:)),'b*', t, y22, 'g*')
xlabel('Time')
ylabel('y_{22}')
legend('Model6','Data')
grid on
axis([-0.2 2 -0.1 0.1])

%when ns = 7
ns = 7;
[U, S, V] = svd(H);
U1 = U([1:m*n],[1:ns]);
S1 = S([1:ns],[1:ns]);
V1 = V([1:n*q],[1:ns]);

On7 = U1 * sqrt(S1);
Cn7 = sqrt(S1) * V1';

On7inv = inv(sqrt(S1)) * U1';
Cn7inv = V1 * inv(sqrt(S1));

C7 = On7([1:m],:);
B7 = Cn7(:,[1:q]);
A7 = On7inv * Htil * Cn7inv;

mev7 = max(abs(eig(A7)));
fprintf('When ns = 7, dominant eval is: %f \n',mev7);

%Plotting impulse response for model7
h = zeros(m,q,n);
x = B7;
for k = 1:n
    h(:,:,k) = C7 * x;
    x = A7 * x;
end

tsim = [1 : n] * ts;
figure()
subplot(211)
plot(tsim,squeeze(h(1,1,:)),'b*', t, y11, 'g*')
xlabel('Time')
ylabel('y_{11}')
legend('Model7', 'Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,1,:)),'b*', t, y21, 'g*')
xlabel('Time')
ylabel('y_{21}')
legend('Model7','Data')
grid on
axis([-0.2 2 -0.1 0.1])

figure()
subplot(211)
plot(tsim,squeeze(h(1,2,:)),'b*', t, y12, 'g*')
xlabel('Time')
ylabel('y_{12}')
legend('Model7','Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,2,:)),'b*', t, y22, 'g*')
xlabel('Time')
ylabel('y_{22}')
legend('Model7','Data')
grid on
axis([-0.2 2 -0.1 0.1])

%when ns = 10
ns = 10;
[U, S, V] = svd(H);
U1 = U([1:m*n],[1:ns]);
S1 = S([1:ns],[1:ns]);
V1 = V([1:n*q],[1:ns]);

On10 = U1 * sqrt(S1);
Cn10 = sqrt(S1) * V1';

On10inv = inv(sqrt(S1)) * U1';
Cn10inv = V1 * inv(sqrt(S1));

C10 = On10([1:m],:);
B10 = Cn10(:,[1:q]);
A10 = On10inv * Htil * Cn10inv;

mev10 = max(abs(eig(A10)));
fprintf('When ns = 10, dominant eval is: %f \n',mev10);

h = zeros(m,q,n);
x = B10;
for k = 1:n
    h(:,:,k) = C10 * x;
    x = A10 * x;
end

tsim = [1 : n] * ts;
figure()
subplot(211)
plot(tsim,squeeze(h(1,1,:)),'b*', t, y11, 'g*')
xlabel('Time')
ylabel('y_{11}')
legend('Model10', 'Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,1,:)),'b*', t, y21, 'g*')
xlabel('Time')
ylabel('y_{21}')
legend('Model10','Data')
grid on
axis([-0.2 2 -0.1 0.1])

figure()
subplot(211)
plot(tsim,squeeze(h(1,2,:)),'b*', t, y12, 'g*')
xlabel('Time')
ylabel('y_{12}')
legend('Model10','Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,2,:)),'b*', t, y22, 'g*')
xlabel('Time')
ylabel('y_{22}')
legend('Model10','Data')
grid on
axis([-0.2 2 -0.1 0.1])

%when ns = 20
ns = 20;
[U, S, V] = svd(H);
U1 = U([1:m*n],[1:ns]);
S1 = S([1:ns],[1:ns]);
V1 = V([1:n*q],[1:ns]);

On20 = U1 * sqrt(S1);
Cn20 = sqrt(S1) * V1';

On20inv = inv(sqrt(S1)) * U1';
Cn20inv = V1 * inv(sqrt(S1));

C20 = On20([1:m],:);
B20 = Cn20(:,[1:q]);
A20 = On20inv * Htil * Cn20inv;

mev20 = max(abs(eig(A20)));
fprintf('When ns = 20, dominant eval is: %f \n',mev20);

h = zeros(m,q,n);
x = B20;
for k = 1:n
    h(:,:,k) = C20 * x;
    x = A20 * x;
end

tsim = [1 : n] * ts;
figure()
subplot(211)
plot(tsim,squeeze(h(1,1,:)),'b*', t, y11, 'g*')
xlabel('Time')
ylabel('y_{11}')
legend('Model20', 'Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,1,:)),'b*', t, y21, 'g*')
xlabel('Time')
ylabel('y_{21}')
legend('Model20','Data')
grid on
axis([-0.2 2 -0.1 0.1])

figure()
subplot(211)
plot(tsim,squeeze(h(1,2,:)),'b*', t, y12, 'g*')
xlabel('Time')
ylabel('y_{12}')
legend('Model20','Data')
grid on
axis([-0.2 2 -0.1 0.1])

subplot(212)
plot(tsim,squeeze(h(2,2,:)),'b*', t, y22, 'g*')
xlabel('Time')
ylabel('y_{22}')
legend('Model20','Data')
grid on
axis([-0.2 2 -0.1 0.1])

fprintf('Since all dominant evals have a magnitude less than 1, all models are asymptotically stable');

%****************  TASK 1, Question 3 and 4  *********************
%Computing Frequency response
wnyq = 20;
w = [0:wnyq/100:wnyq];%w in Hertz
Hf6 = zeros(m,q,length(w));
Hf7 = zeros(m,q,length(w));
Hf10 = zeros(m,q,length(w));
Hf20 = zeros(m,q,length(w));

for k = 1:length(w)
    Hf6(:,:,k) = C6 * inv(exp(1j*2*pi*w(k)*ts)*eye(6) - A6) * B6;
    Hf7(:,:,k) = C7 * inv(exp(1j*2*pi*w(k)*ts)*eye(7) - A7) * B7;
    Hf10(:,:,k) = C10 * inv(exp(1j*2*pi*w(k)*ts)*eye(10) - A10) * B10;
    Hf20(:,:,k) = C20 * inv(exp(1j*2*pi*w(k)*ts)*eye(20) - A20) * B20;
end

%Computing emperical frequency response
y11f = fft(y11)./fft(u1);
N = length(y11f);
om = [0:N-1]/(ts*N); %%%% frequency vector in hertz
y21f = fft(y21)./fft(u1);
y12f = fft(y12)./fft(u2);
y22f = fft(y22)./fft(u2);

%Plotting abs(Freq response of channel (1,1))
figure()
plot(w, abs(squeeze(Hf6(1,1,:))), 'b-')
hold on
plot(w, abs(squeeze(Hf7(1,1,:))), 'r-')
plot(w, abs(squeeze(Hf10(1,1,:))), 'g-')
plot(w, abs(squeeze(Hf20(1,1,:))), 'm-')
plot(om,abs(y11f), 'k--')
ylabel('|H(w)| (1,1)')
xlabel('Frequency (Hz)')
xlim([0 20])
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response magnitude of (1,1) channel')
grid on
hold off

%Plotting abs(Freq response) of channel (2,1)
figure()
plot(w, abs(squeeze(Hf6(2,1,:))), 'b-')
hold on
plot(w, abs(squeeze(Hf7(2,1,:))), 'r-')
plot(w, abs(squeeze(Hf10(2,1,:))), 'g-')
plot(w, abs(squeeze(Hf20(2,1,:))), 'm-')
plot(om,abs(y21f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('|H(w)| (2,1)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response magnitude of (2,1) channel')
grid on
hold off

%Plotting abs(Freq response) of channel (1,2)
figure()
plot(w, abs(squeeze(Hf6(1,2,:))), 'b-')
hold on
plot(w, abs(squeeze(Hf7(1,2,:))), 'r-')
plot(w, abs(squeeze(Hf10(1,2,:))), 'g-')
plot(w, abs(squeeze(Hf20(1,2,:))), 'm-')
plot(om,abs(y12f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('|H(w)| (1,2)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response magnitude of (1,2) channel')
grid on
hold off

%Plotting abs(Freq response) of channel (2,2)
figure()
plot(w, abs(squeeze(Hf6(2,2,:))), 'b-')
hold on
plot(w, abs(squeeze(Hf7(2,2,:))), 'r-')
plot(w, abs(squeeze(Hf10(2,2,:))), 'g-')
plot(w, abs(squeeze(Hf20(2,2,:))), 'm-')
plot(om,abs(y22f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('|H(w)| (2,2)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response magnitude of (2,2) channel')
grid on
hold off

%Plotting angle(Freq response of channel (1,1))
figure()
plot(w, angle(squeeze(Hf6(1,1,:))), 'b-')
hold on
plot(w, angle(squeeze(Hf7(1,1,:))), 'r-')
plot(w, angle(squeeze(Hf10(1,1,:))), 'g-')
plot(w, angle(squeeze(Hf20(1,1,:))), 'm-')
plot(om,angle(y11f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('angle(H(w)) (1,1)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response angle of (1,1) channel')
grid on
hold off

%Plotting angle(Freq response) of channel (2,1)
figure()
plot(w, angle(squeeze(Hf6(2,1,:))), 'b-')
hold on
plot(w, angle(squeeze(Hf7(2,1,:))), 'r-')
plot(w, angle(squeeze(Hf10(2,1,:))), 'g-')
plot(w, angle(squeeze(Hf20(2,1,:))), 'm-')
plot(om,angle(y21f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('angle(H(w))(2,1)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response angle of (2,1) channel')
grid on
hold off

%Plotting angle(Freq response) of channel (1,2)
figure()
plot(w, angle(squeeze(Hf6(1,2,:))), 'b-')
hold on
plot(w, angle(squeeze(Hf7(1,2,:))), 'r-')
plot(w, angle(squeeze(Hf10(1,2,:))), 'g-')
plot(w, angle(squeeze(Hf20(1,2,:))), 'm-')
plot(om,angle(y12f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('angle(H(w)) (1,2)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response angle of (1,2) channel')
grid on
hold off

%Plotting angle(Freq response) of channel (2,2)
figure()
plot(w, angle(squeeze(Hf6(2,2,:))), 'b-')
hold on
plot(w, angle(squeeze(Hf7(2,2,:))), 'r-')
plot(w, angle(squeeze(Hf10(2,2,:))), 'g-')
plot(w, angle(squeeze(Hf20(2,2,:))), 'm-')
plot(om,angle(y22f), 'k--')
xlabel('Frequency (Hz)')
xlim([0 20])
ylabel('angle(H(w)) (2,2)')
legend('Model6','Model7','Model10','Model20','Data')
title('Freq Response phase of (2,2) channel')
grid on
hold off



