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

ev7 = eig(A7);

up = sum(diag(S1));
lo = 0;
tol = 1e-6;

func = @Hinf_cont;
[H_val, f_val] = Hinf_dis(func, A7, B7, C7, zeros(2), up, lo, tol, ts);
%H_val is the H_inf norm (0.4693 for our ns = 7 model), and it occurs at a
%frequency f_val (71.1487 rad/s).