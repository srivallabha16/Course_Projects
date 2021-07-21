function [Hinf, freq] = Hinf_dis(Hcon,A, B, C, D, up, lo, tol, tsam)
%INPUTS: Function to compute H_inf for continuous time systems,
%        discrete time matrices - A, B, C, D
%        upper limit, lower limit, and tolerance for gamma
%RETURNS the H_inf norm and the discrete time frequency at which it ocuurs.
Ac = -inv(eye(size(A)) + A) * (eye(size(A)) - A);
Bc = sqrt(2) * inv(eye(size(A)) + A) * B;
Cc = sqrt(2) * C * inv(eye(size(A)) + A);
Dc = D - C * inv(eye(size(A)) + A) * B;

[Hinf, freq] = Hcon(Ac, Bc, Cc, Dc, up, lo, tol);
freq = 1/tsam * angle((1 + 1i * freq)/ (1 - 1i * freq));   

end

