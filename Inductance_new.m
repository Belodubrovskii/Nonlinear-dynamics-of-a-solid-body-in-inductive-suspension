clear; clc; close all
% 

a = 2;  % radius of bottom ring
b = 2;   
c = 2;  % distance

global M g Self_inductance Im Omega Mu_fric Resistance %dL12dZ L12

Mu = 12.57*10^(-7);  % permeability of free space ÑÈ
Mu0 = 1;  % ÑÃÑ
Mu1 = 1; %conductor permeability
R0 = 0.006;  % ring radius
R = 0.006;   % coil radius
Im = 0.2;
Omega = 4*10^3*2*pi;
r0 = 1.5e-5 ; % ring section radius   
g = 9.8;
Mu_fric = 0.0001 ;
S_ring = pi*r0^2;
l_ring = 2*pi*R0;
Resistance = 0.0171*S_ring/l_ring;
M = S_ring*l_ring*8960;

syms alpha0 beta0 ez er tetta X Y Z
assume(alpha0,'real'); assume(beta0,'real'); assume(tetta,'real');
assume(X,'real'); assume(Y,'real'); assume(Z,'real');
% !(er, e_phi, ez)
n = [cos(alpha0); sin(alpha0)*cos(beta0); sin(alpha0)*sin(beta0)];
ez = [0; 0; 1];

Qx = Turn([1; 0; 0],alpha0);
Qz = Turn([0; 0; 1],beta0);
Q = Qz*Qx;

b1 = Q(:,1);
b2 = Q(:,2);

e = b2*cos(tetta) - b1*sin(tetta);

i = [1; 0; 0]; j = [0; 1; 0]; %k = [0; 0; 1];
r = [X; Y; Z] + R0*(b1*cos(tetta) + b2*sin(tetta));  %vector to point on the ring
r_projection = [r(1); r(2); 0];

cos_phi = r(1)/sqrt(r_projection.'*r_projection);
sin_phi = r(2)/sqrt(r_projection.'*r_projection);
e_phi = j*cos_phi - i*sin_phi;
% disp(e_phi);

r_projection_length = norm(r_projection);
z = r(3);

k_square = 4*R*r_projection_length /((r_projection_length + R)^2 + z^2);
[K,E] = ellipke(k_square);

Vect_potential = Mu / (pi*sqrt(k_square))*(R/r_projection_length)^(1/2) * ((1 - k_square/2)*K - E);
funSym = Vect_potential*R0*dot(e_phi,e);
funSym = simplify(funSym,'ignoreAnalyticConstraints',true);

if ~exist('L12.m','file')
    matlabFunction(funSym,'file','L12.m','vars',[alpha0 beta0 tetta X Y Z]);
    error('Please change ellipticK/E functions calls with ellipke')
end

dfunSymdZ = simplify(diff(funSym,Z),'ignoreAnalyticConstraints',true);
dfunSymdY = diff(funSym,Y);
dfunSymdX = diff(funSym,X);
dfunSymdalpha0 = diff(funSym,alpha0);
dfunSymdbeta0 = diff(funSym,beta0);
if ~exist('dL12dZ.m','file')
    matlabFunction(dfunSymdZ,'file','dL12dZ.m','vars',[alpha0 beta0 tetta X Y Z]);
    error('Please change ellipticK/E functions calls with ellipke')
end

dL12dY = matlabFunction(dfunSymdY,'vars',[alpha0 beta0 tetta X Y Z]);
dL12dX = matlabFunction(dfunSymdX,'vars',[alpha0 beta0 tetta X Y Z]);
dL12dalpha0 = matlabFunction(dfunSymdalpha0,'vars',[alpha0 beta0 tetta X Y Z]);
dL12dbeta0 = matlabFunction(dfunSymdbeta0,'vars',[alpha0 beta0 tetta X Y Z]);

alpha0Num = 0; beta0Num = 0; XNum = 0; YNum = 0; ZNum = 2.5e-5;
Mutual_inductance = integral(@(x) L12(alpha0Num,beta0Num,x,XNum,YNum,ZNum), 0, 2*pi);


%Self_inductance = pi*Mu1*R + Mu0*4*pi*R*(log10(8*R/r0) - 2);
Self_inductance = Mu*R0*(log(8*R0/r0) - 2 + 1/4);



% disp(Mutual_inductance);



%==========================================================================
%teoretical expression
k = 4*a*b/((a + b)^2 + c^2);
[K1,E1] = ellipke(k);
Mutual_inductance_th = 2*Mu/sqrt(k)*(a*b)^(1/2)*((1-k/2)*K1-E1);

% disp(Mutual_inductance_th);
%==========================================================================

% disp(Self_inductance);

%============ ODE SOLVER ===================================
%options = odeset('OutputFcn',@odeplot,'OutputSel',1);
[t,z] = ode45(@(t,z) odeFun(t,z),[0 1],[2.5e-5, 0, 0]);
figure(); plot(t,z(:,1));
figure(); plot(t,z(:,3));
figure(); plot(t,z(:,2));

%figure(); plot(t,dL12dZ(0,0,0,0,0,z(:,1)).*z(:,2)*Im.*cos(Omega*t));
%figure(); plot(t,1/2*Omega*L12(0,0,0,0,0,z(:,1))*Im.*sin(Omega*t));


%============================================================
function dydt = odeFun(t,z)
global M g Self_inductance Im Omega Mu_fric  Resistance %dL12dZ L12
alpha0Num = 0; beta0Num = 0; XNum = 0; YNum = 0;
dydt = zeros(3,1);
dydt(1) = z(2);
% dL12 = integral(@(x) dL12dZ(alpha0Num,beta0Num,x,XNum,YNum,z(1)), 0, 2*pi);
% L12_fun = integral(@(x) L12(alpha0Num,beta0Num,x,XNum,YNum,z(1)), 0, 2*pi);
dL12 = 2*pi*dL12dZ(alpha0Num,beta0Num,0,XNum,YNum,z(1));
L12_fun = 2*pi*L12(alpha0Num,beta0Num,0,XNum,YNum,z(1));
dydt(2) = 1/2/M*dL12*z(3)*Im*cos(Omega*t) - Mu_fric*z(2)/M - g ;
dydt(3) = (-1/2*dL12*z(2)*Im*cos(Omega*t) + 1/2*Omega*L12_fun*Im*sin(Omega*t) - Resistance*z(3)) / Self_inductance;

end

function Q = Turn(vector, angle)
mxE = zeros(3); ort = [1; 0; 0];
for i = 1:3
    mxE(:,i) = cross(vector,circshift(ort,i-1));
end

Q = (1 - cos(angle))*kron(vector.', vector) + cos(angle)*eye(3) + sin(angle)*mxE;
end



    