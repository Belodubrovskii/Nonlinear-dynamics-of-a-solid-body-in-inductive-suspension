clear; clc; close all

Mu = 12.57*10^(-7);  % permeability of free space СИ
% R0 = 0.0015;  % ring radius
% R = 0.002;   % coil radius
Im = 0.2;
Omega = 4*10^3*2*pi;
% r0 = 1e-4 ; % ring section radius   
g = 9.8;
Mu_fric = 0.0001 ;

syms alpha0 beta0 ez er tetta X Y Z R0 R r0
assume(alpha0,'real'); assume(beta0,'real'); assume(tetta,'real');
assume(X,'real'); assume(Y,'real'); assume(Z,'real');

S_ring = pi*r0^2;
l_ring = 2*pi*R0;
S_ring_mm = S_ring*10^6;
Resistance = 0.0171*l_ring/S_ring_mm;
M = S_ring*l_ring*8960;

MFun = matlabFunction(M,'vars',[r0 R0]);

n = [cos(alpha0); sin(alpha0)*cos(beta0); sin(alpha0)*sin(beta0)];
ez = [0; 0; 1];

Qx = Turn([1; 0; 0],alpha0);
Qz = Turn([0; 0; 1],beta0);
Q = Qz*Qx;

b1 = Q(:,1);
b2 = Q(:,2);

pretty(b1)
pretty(b2)

e = b2*cos(tetta) - b1*sin(tetta);

i = [1; 0; 0]; j = [0; 1; 0]; %k = [0; 0; 1];
r = [X; Y; Z] + R0*(b1*cos(tetta) + b2*sin(tetta));  %vector to point on the ring
r_projection = [r(1); r(2); 0];

cos_phi = r(1)/sqrt(r_projection.'*r_projection);
sin_phi = r(2)/sqrt(r_projection.'*r_projection);
e_phi = j*cos_phi - i*sin_phi;

r_projection_length = norm(r_projection);
z = r(3);

k_square = 4*R*r_projection_length /((r_projection_length + R)^2 + z^2);
[K,E] = ellipke(k_square);

Vect_potential = Mu / (pi*sqrt(k_square))*(R/r_projection_length)^(1/2) * ((1 - k_square/2)*K - E);
funSym = Vect_potential*R0*2*pi*dot(e_phi,e);
funSym = simplify(funSym,'ignoreAnalyticConstraints',true);

dfunSymdZ = simplify(diff(funSym,Z),'ignoreAnalyticConstraints',true);


L12Fun = matlabFunction(funSym,'vars',[alpha0 beta0 tetta r0 R0 R X Y Z]);
dL12dZFun =  matlabFunction(dfunSymdZ,'vars',[alpha0 beta0 tetta r0 R0 R X Y Z]);

alpha0Num = 0; beta0Num = 0; TettaNum = 0; XNum = 0; YNum = 0; ZNum = 2.5e-5;

L12Fun = @(r0,R0,R) L12Fun(alpha0Num,beta0Num, TettaNum,r0,R0,R,XNum,YNum,ZNum);
dL12dZFun = @(r0,R0,R) dL12dZFun(alpha0Num,beta0Num, TettaNum,r0,R0,R,XNum,YNum,ZNum);


L22 = Mu*R0*(log(8*R0/r0) - 2 + 1/4);
L22Fun = matlabFunction(L22,'vars',[r0 R0]);

ResistanceFun = matlabFunction(Resistance,'vars',[r0 R0]);

A = @(r0,R0,R) 1/2*Omega*Im.*L12Fun(r0,R0,R)/sqrt((L22Fun(r0,R0)*Omega).^2 + (ResistanceFun(r0,R0)).^2);

% Complex_amplitude = s + d*i
s = @(r0,R0,R) -1/2*L12Fun(r0,R0,R)*Im*Omega^2*L22Fun(r0,R0)/((L22Fun(r0,R0)*Omega).^2 ...
    + (ResistanceFun(r0,R0)).^2);
d = @(r0,R0,R) -1/2*L12Fun(r0,R0,R)*ResistanceFun(r0,R0)*Im*Omega...
    /((L22Fun(r0,R0)*Omega).^2 + (ResistanceFun(r0,R0)).^2);

Phi = @(r0,R0,R) atan2(d(r0,R0,R),s(r0,R0,R));

% L22Fun(r0,R0)*Omega./ResistanceFun(r0,R0)

Func = @(Z) (MFun(Z(1),Z(2))*g) ./ ...
    (1/2*Im*A(Z(1), Z(2), Z(3)).*cos(Phi(Z(1),Z(2),Z(3))).*dL12dZFun(Z(1),Z(2),Z(3)));

% Func2 = @(Z) (1/2*Im*A(Z(1), Z(2), Z(3)).*cos(Phi(Z(1), Z(2))).*dL12dZFun(Z(1),Z(2),Z(3)))...
% / (MFun(Z(1),Z(2))*g) - 1;

%=============================================================================
%Подбор параметров r0 R0 R

x0 = [1e-5, 0.002, 0.003];
lb = [1e-6, 0.0005, 0.0005];
ub = [1e-4, 0.006, 0.006];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x,fval] = fmincon(Func,x0,[],[],[],[], lb, ub,[],options);
y = vpa(x,7);



% disp([x,fval]);


% tettaNum = 0;
% zz = linspace(1e-7,1e-2,15);
% L12_fun = 2*pi*L12Fun(alpha0Num,beta0Num,tettaNum,XNum,YNum,zz);
% A_fun = 1/2*Omega*Im*L12_fun/sqrt((L22*Omega)^2 + Resistance^2);
% plot(zz,A_fun);




function Q = Turn(vector, angle)
mxE = zeros(3); ort = [1; 0; 0];
for i = 1:3
    mxE(:,i) = cross(vector,circshift(ort,i-1));
end

Q = (1 - cos(angle))*kron(vector.', vector) + cos(angle)*eye(3) + sin(angle)*mxE;
end
