%********************   TASK 2  *************************
clear
close all
clc

load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[~, mi] = max(u1>0); %%% find index where pulse occurs

load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11(1:mi-1));
y12 = y12 - mean(y12(1:mi-1));
y21 = y21 - mean(y21(1:mi-1));
y22 = y22 - mean(y22(1:mi-1));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);

%****************  TASK 2, Question 1 and 2  *********************

ts = 1/40; %%%% sample period

m = 2; %number of output channels
q = 2; %number of input channels
n = 100; % n s.t size(H) = mn X qn

ind_off = 41;
H = zeros(m*n, q*n);
Htil = zeros(m*n, q*n);

%To create H and H tilda matrices
r = 1;
c = 1;
for k = 1 : n
    ind = ind_off + k + (0:n-1);     
    
    H(r,c:q:c+q*(n-1)) = y11(ind);
    H(r+1,c:q:c+q*(n-1)) = y21(ind);
    H(r,c+1:q:c+1+q*(n-1)) = y12(ind);
    H(r+1,c+1:q:c+1+q*(n-1)) = y22(ind);
    
    Htil(r,c:q:c+q*(n-1)) = y11(ind+1);
    Htil(r+1,c:q:c+q*(n-1)) = y21(ind+1);
    Htil(r,c+1:q:c+1+q*(n-1)) = y12(ind+1);
    Htil(r+1,c+1:q:c+1+q*(n-1)) = y22(ind+1);
    
    r = r + m;
end

%when ns = 7
ns = 7;
[U, S, V] = svd(H);
U1 = U(1:m*n, 1:ns);
S1 = S(1:ns, 1:ns);
V1 = V(1:n*q, 1:ns);

On7 = U1 * sqrt(S1);
Cn7 = sqrt(S1) * V1';

On7inv = inv(sqrt(S1)) * U1';
Cn7inv = V1 * inv(sqrt(S1));

C7 = On7([1:m],:);
B7 = Cn7(:,[1:q]);
A7 = On7inv * Htil * Cn7inv;

M = [A7 B7; -C7 zeros(2)];
N = [eye(7) zeros(7,2); zeros(2,7) zeros(2)];%disp(diag(D))%Resulted in D(1) and D(8) as inf
[~, D] = eig(M, N);


%2,3,4,5,6 elements of the diagonal are not inf and are in descending
%order. Elements 4 and 5 of the diagonal are equal
z1 = D(2,2); z2 = D(3,3); z3 = D(4,4); z4 = D(5,5); z5 = D(6,6);
disp('The finite transmission zeros are\n')
fprintf('The five finite transmission zeros are:\n')
fprintf('z1 is %f + %fi and |z1| is %f \n',real(z1),imag(z1),abs(z1))
fprintf('z2 is %f + %fi and |z2| is %f \n',real(z2),imag(z2),abs(z2))
fprintf('z3 is %f + %fi and |z3| is %f \n',real(z3),imag(z3),abs(z3))
fprintf('z4 is %f + %fi and |z4| is %f \n',real(z4),imag(z4),abs(z4))
fprintf('z5 is %f + %fi and |z5| is %f \n',real(z5),imag(z5),abs(z5))

eigA7 = eig(A7);
ut = 0:2*pi/100:2*pi;
figure()
plot(eigA7,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
title('Poles and Zeros')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-3 2.5 -1 1])
axis equal
grid on
hold off

%****************  TASK 2, Question 3  *********************
%Converting discrete evals to continuous evals
eigc_real = 1/ts * log(abs(eigA7));
eigc_im = 1/ts * angle(eigA7);
eigc = eigc_real + 1j * eigc_im;

fprintf('\n lambda_continuous \n');
disp(eigc);

disp('eigc([1,2,3,4]) have Re(.) < 0')
disp('There are two complex conjugate pairs with negative real parts (2 damped oscillators)')
fprintf('\n Oscillator 1 has a freq %f rad/sec = %f hertz\n',imag(eigc(1)),imag(eigc(1))/(2*pi));

%****************  TASK 2, Question 4  *********************

%For 1-1 channel the system is 
B11 = B7(:,1);
C11 = C7(1,:);
M = [A7 B11; -C11 zeros(1)];
N = [eye(7) zeros(7,1);zeros(1,7) zeros(1)];
[~, D] = eig(M,N);
z11 = diag(D);
z11 = z11([2,3,4,5,6,7]);
disp('Finite zeros for 1-1 channel:')
disp(z11)
lic = num2cell(z11);
[~, z2, z3, z4, z5, z6]=deal(lic{:});
figure()
plot(eigA7,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
title('Poles and Zeros 1-1 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 1-2 channel the system is 
B12 = B7(:,2);
C12 = C7(1,:);
M = [A7 B12; -C12 zeros(1)];
N = [eye(7) zeros(7,1);zeros(1,7) zeros(1)];
[~, D] = eig(M,N);
z12 = diag(D);
z12= z12([2,3,4,5,6,7]);
disp('Finite zeros for 1-2 channel:')
disp(z12)
lic = num2cell(z12);
[~, z2, z3, z4, z5, z6]=deal(lic{:});
figure()
plot(eigA7,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
title('Poles and Zeros 1-2 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 2-1 channel
B21 = B7(:,1);
C21 = C7(2,:);
M = [A7 B21; -C21 zeros(1)];
N = [eye(7) zeros(7,1);zeros(1,7) zeros(1)];
[~, D] = eig(M,N);
z21 = diag(D);
z21 = z21([2,3,4,5,6,7]);
disp('Finite zeros for 2-1 channel:')
disp(z21)
lic = num2cell(z21);
[z1, z2, z3, z4, z5, z6]=deal(lic{:});
figure()
plot(eigA7,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
title('Poles and Zeros 2-1 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 2-2 channel
B22 = B7(:,2);
C22 = C7(2,:);
M = [A7 B22; -C22 zeros(1)];
N = [eye(7) zeros(7,1);zeros(1,7) zeros(1)];
[~, D] = eig(M,N);
z22 = diag(D);
z22 = z22([2,3,4,5,6,7]);
disp('Finite zeros for 2-2 channel:')
disp(z22)
lic = num2cell(z22);
[~, z2, z3, z4, z5, z6]=deal(lic{:});
figure()
plot(eigA7,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
title('Poles and Zeros 2-2 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%Hankel Matrix for 1-1 channel
m = 1; q =1;
H11 = zeros(m*n, q*n);
H12 = zeros(m*n, q*n);
H21 = zeros(m*n, q*n);
H22 = zeros(m*n, q*n);
r = 1;
c = 1;
for k = 1 : n
    ind = ind_off + k + [0:n-1]; 
    
    H11(r,:) = y11(ind);
    H12(r,:) = y12(ind);
    H21(r,:) = y21(ind);
    H22(r,:) = y22(ind);

    r = r + m;
end

figure()
subplot(411)
temp = svd(H11);
disp('The first 5 Hankel singular values for the 1-1 channel are')
disp(temp([1:5]));
plot(temp,'bo','DisplayName','H11');
ylabel('H11 sing vals');
xlabel('Index');xlim([0,20]);
legend
hold on

subplot(412)
temp = svd(H12);
disp('The first 5 Hankel singular values for the 1-2 channel are')
disp(temp([1:5]));
plot(temp,'m+','DisplayName','H12');
ylabel('H12 sing vals');
xlabel('Index');xlim([0,20]);
legend
hold on

subplot(413)
temp = svd(H21);
disp('The first 5 Hankel singular values for the 2-1 channel are')
disp(temp([1:5]));
plot(temp,'k*','DisplayName','H21');
ylabel('H21 sing vals');
xlabel('Index');
xlim([0,20]);
legend
hold on

subplot(414)
temp = svd(H22);
disp('The first 5 Hankel singular values for the 1-2 channel are')
disp(temp([1:5]));
plot(temp,'rd','DisplayName','H22');
ylabel('H22 sing vals');
xlabel('Index');
xlim([0,20]);
legend
hold on

%****************  TASK 2, Question 5  *********************

%Generating A,B, and C for ns = 8
m = 2; %number of output channels
q = 2; %number of input channels
ns = 8;
[U, S, V] = svd(H);
U1 = U([1:m*n],[1:ns]);
S1 = S([1:ns],[1:ns]);
V1 = V([1:n*q],[1:ns]);

On8 = U1 * sqrt(S1);
Cn8 = sqrt(S1) * V1';

On8inv = inv(sqrt(S1)) * U1';
Cn8inv = V1 * inv(sqrt(S1));

C8 = On8([1:m],:);
B8 = Cn8(:,[1:q]);
A8 = On8inv * Htil * Cn8inv;
eigA8 = eig(A8);

%For 1-1 channel the system is 
B11 = B8(:,1);
C11 = C8(1,:);
M = [A8 B11; -C11 zeros(1)];
N = [eye(8) zeros(8,1);zeros(1,8) zeros(1)];
[~, D] = eig(M,N);
z11 = diag(D);
z11 = z11([2,3,4,5,6,7,8]);
disp('Finite zeros for 1-1 channel:')
disp(z11)
lic = num2cell(z11);
[z1 z2 z3 z4 z5 z6 z7]=deal(lic{:});
figure()
plot(eigA8,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
%plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
plot(real(z7),imag(z7),'yo','DisplayName',['|z7| is ',num2str(abs(z7))]);
title('ns = 8, Poles and Zeros 1-1 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 1-2 channel the system is 
B12 = B8(:,1);
C12 = C8(2,:);
M = [A8 B12; -C12 zeros(1)];
N = [eye(8) zeros(8,1);zeros(1,8) zeros(1)];
[V,D] = eig(M,N);
%disp(diag(D))%Resulted in D(1) and D(9) as inf
z12 = diag(D);
z12 = z12([2,3,4,5,6,7,8]);
disp('Finite zeros for 1-2 channel:')
disp(z12)
lic = num2cell(z12);
[z1 z2 z3 z4 z5 z6 z7]=deal(lic{:});
figure()
plot(eigA8,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
%plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
plot(real(z7),imag(z7),'yo','DisplayName',['|z7| is ',num2str(abs(z7))]);
title('ns = 8, Poles and Zeros 1-2 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 2-1 channel the system is 
B21 = B8(:,2);
C21 = C8(1,:);
M = [A8 B21; -C21 zeros(1)];
N = [eye(8) zeros(8,1);zeros(1,8) zeros(1)];
[V,D] = eig(M,N);
%disp(diag(D))%Resulted in D(1) and D(9) as inf
z21 = diag(D);
z21 = z21([2,3,4,5,6,7,8]);
disp('Finite zeros for 2-1 channel:')
disp(z21)
lic = num2cell(z21);
[z1 z2 z3 z4 z5 z6 z7]=deal(lic{:});
figure()
plot(eigA8,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
%plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
plot(real(z7),imag(z7),'yo','DisplayName',['|z7| is ',num2str(abs(z7))]);
title('ns = 8, Poles and Zeros 2-1 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off

%For 2-2 channel the system is 
B22 = B8(:,2);
C22 = C8(2,:);
M = [A8 B22; -C22 zeros(1)];
N = [eye(8) zeros(8,1);zeros(1,8) zeros(1)];
[V,D] = eig(M,N);
%disp(diag(D))%Resulted in D(1) and D(9) as inf
z22 = diag(D);
z22 = z22([2,3,4,5,6,7,8]);
disp('Finite zeros for 2-2 channel:')
disp(z22)
lic = num2cell(z22);
[z1 z2 z3 z4 z5 z6 z7]=deal(lic{:});
figure()
plot(eigA8,'x','DisplayName','Poles');
hold on
plot(cos(ut),sin(ut),'k-','DisplayName','Unit Circle')
%plot(real(z1),imag(z1),'bo','DisplayName',['|z1| is ',num2str(abs(z1))]);
plot(real(z2),imag(z2),'go','DisplayName',['|z2| is ',num2str(abs(z2))]);
plot(real(z3),imag(z3),'ko','DisplayName',['|z3| is ',num2str(abs(z3))]);
plot(real(z4),imag(z4),'mo','DisplayName',['|z4| is ',num2str(abs(z4))]);
plot(real(z5),imag(z5),'co','DisplayName',['|z5| is ',num2str(abs(z5))]);
plot(real(z6),imag(z6),'ro','DisplayName',['|z6| is ',num2str(abs(z6))]);
plot(real(z7),imag(z7),'yo','DisplayName',['|z7| is ',num2str(abs(z7))]);
title('ns = 8, Poles and Zeros 2-2 channel')
legend
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-1 1 -1 1])
axis equal
grid on
hold off
