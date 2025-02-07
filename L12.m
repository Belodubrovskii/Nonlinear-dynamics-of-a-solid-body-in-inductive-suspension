function funSym = L12(alpha0,beta0,tetta,X,Y,Z)
%L12
%    FUNSYM = L12(ALPHA0,BETA0,TETTA,X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    12-Apr-2020 15:12:49

t3 = sin(beta0);
t4 = cos(tetta);
t5 = cos(alpha0);
t6 = cos(beta0);
t7 = sin(tetta);
t9 = X.*5.0e2;
t10 = t4.*t6.*3.0;
t11 = t3.*t5.*t7.*3.0;
t2 = t9+t10-t11;
t13 = Y.*5.0e2;
t14 = t3.*t4.*3.0;
t15 = t5.*t6.*t7.*3.0;
t8 = t13+t14+t15;
t12 = t2.^2;
t16 = t8.^2;
t17 = t12+t16;
t20 = sqrt(t17);
t18 = t20+3.0;
t22 = Z.*5.0e2;
t23 = sin(alpha0);
t24 = t7.*t23.*3.0;
t19 = t22+t24;
t21 = t18.^2;
t25 = t19.^2;
t26 = t21+t25;
t27 = 1.0./t26;
t28 = t20.*t27.*1.2e1;
[t29, t30] = ellipke(t28);
funSym = (sqrt(t26).*(t30 + t29.*(t20.*t27.*6.0-1.0)).*(t5.*3.0-X.*t3.*t7.*5.0e2+Y.*t6.*t7.*5.0e2+X.*t4.*t5.*t6.*5.0e2+Y.*t3.*t4.*t5.*5.0e2).*(-3.771e-9))./(t17.*pi);
