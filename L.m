function L
%           Построение ДА РЛС
Rmax = 9.5783e+04;          %максимальная дальность обнаружения с учетом затухания
%угды для апроксимации ДН
z0 = pi/180*1.8;       %направление излучения
zmin = pi/180 * 0.1;   %минимальный угол излучения
zmax = pi/180 * 5;     %максимальный угол излучния

ha = 16;                %высота антены
hl=1.05;               %длина волны
k0=2*pi/hl;       %волновое число 
E0=1/(36*pi) * 10^(-9); %электрическая постоянная вакуума
E2=6;    %диэлектрическая проницаемость подстилающей поверхности
E1=1;    %диэлектрическая проницаемость атмосферы (так во всех вар-х)
u0=4*pi*10^(-7);        %магнитная постоянная вакуума
u1=1;      %магнитные пронецаемости (так во всех вар-х)
u2=1; 
q1=0;       %удельная проводимость среды (у всех так)
q2=1.5;     %удельная проводимость среды (по вар-м)

%разбивка на дискреты углов (от 0 до 90 град)
delta=0.00001;
b=pi/180*90; 
z=0:delta:b; 

for i=1:numel(z)
if z(i)>zmax
F(i)=0;
end
if z(i)>=z0 & z(i)<zmax
F(i)=sin(z0)*csc(z(i));
end
if z(i)>=zmin & z(i)<z0
F(i)=1;
end
if z(i)<zmin 
F(i)=1-(zmin-z(i))/zmin;
end
end
figure
polar(z,F)
for i=1:numel(z)
X(i)=cos(z(i))*F(i);
Y(i)=sin(z(i))*F(i);
end
figure
plot(X,Y)

for i=1:numel(z)
X(i)=cos(z(i))*F(i)*Rmax;
Y(i)=sin(z(i))*F(i)*Rmax;
end
figure
plot(X,Y)

w=(2*pi*3*10^8)/hl;   %круговая частота

Ea1=E0*E1-(q1/w)*sqrt(-1);
Ea2=E0*E2-(q2/w)*sqrt(-1);

ua1 = u0*u1;
ua2 = u0*u2;

k1=w*sqrt(Ea1*ua1);
k2=w*sqrt(Ea2*ua2);

for j=1:numel(z)
Rv = (Ea2*k1*sin(z(j))-Ea1*sqrt(k2^2-k1^2*(cos(z(j)))^2))/(Ea2*k1*sin(z(j))+Ea1*sqrt(k2^2-k1^2*(cos(z(j)))^2));
Bv = angle(Rv);
L(j)=sqrt(1 + (abs(Rv))^2 + 2*abs(Rv)*cos(2*k0*ha*sin(z(j)) + Bv));
end
figure
plot(z,L)

for i=1:numel(z)
    r(i) = Rmax * F(i) * L(i);
    X_(i) = r(i) * cos(z(i));
    Y_(i) = r(i) * sin(z(i));
end
figure
plot(X,Y,X_,Y_)

dndz = -6*10^-8; %градиент коэффициента преломления радиоволны в тропосфере
Ro = 6370 * 10^3;   %радиус земли

Ria = Ro / (1 + Ro * dndz)    %эквивалентный радиус земли

delta_h = 1000;  %разница высот для сетки

%
Rz = Ro;
Rz2 = Ria;
dh = delta_h;
for i=1:numel(z)
X1(i)=cos(z(i))*Rmax * L(i)*F(i);
Y1(i)=ha + sin(z(i))*F(i)*Rmax * L(i); 
X111(i)=cos(z(i)+pi/180*0)*(Rmax*1) * L(i) * F(i);
Y111(i)=ha + sin(z(i)+pi/180*0)*F(i)*(Rmax*1) * L(i); 
Xm2(i) = i/(numel(z)) * 10*10^4;
Y2(i)=-Xm2(i)^2/(2*Rz);
Xm22(i) = i/(numel(z)) * 1.1*10^5;
Y22(i)=-Xm22(i)^2/(2*Rz2);
Xm3(i) = i/(numel(z)) * 1.13*10^5;
Y3(i)=dh*0.325-Xm3(i)^2/(2*Rz);
Y33(i)=dh*0.325-Xm3(i)^2/(2*Rz2);
Xm30(i) = i/(numel(z)) * 1.13*10^5;
Y30(i)=dh-Xm3(i)^2/(2*Rz);
Y330(i)=dh-Xm3(i)^2/(2*Rz2);
Xm4(i) = i/(numel(z)) * 1.6*10^5;
Y4(i)=2*dh-Xm4(i)^2/(2*Rz);
Y44(i)=2*dh-Xm4(i)^2/(2*Rz2);
Xm5(i) = i/(numel(z)) * 1.95*10^5;
Y5(i)=3*dh-Xm5(i)^2/(2*Rz);
Y55(i)=3*dh-Xm5(i)^2/(2*Rz2);
Xm6(i) = i/(numel(z)) * 2.25*10^5;
Y6(i)=4*dh-Xm6(i)^2/(2*Rz);
Y66(i)=4*dh-Xm6(i)^2/(2*Rz2);
Xm7(i) = i/(numel(z)) * 2.49*10^5;
Y7(i)=5*dh-Xm7(i)^2/(2*Rz);
Y77(i)=5*dh-Xm7(i)^2/(2*Rz2);
Xm8(i) = i/(numel(z)) * 2.49*10^5;
Y8(i)=6*dh-Xm8(i)^2/(2*Rz);
Y88(i)=6*dh-Xm8(i)^2/(2*Rz2);
Xm9(i) = i/(numel(z)) * 2.5*10^5;
Y9(i)=7*dh-Xm9(i)^2/(2*Rz);
Y99(i)=7*dh-Xm9(i)^2/(2*Rz2);
Xm10(i)= i/(numel(z)) * 2.2*10^5;
Y10(i)=ha+tan(pi/180*1)*Xm10(i);
Xm11(i)= i/(numel(z)) * 1.55*10^5;
Y11(i)=ha+tan(pi/180*2)*Xm11(i);
Xm12(i)= i/(numel(z)) * 0.95*10^5;
Y12(i)=ha+tan(pi/180*4)*Xm12(i);
Xm13(i)= i/(numel(z)) * 0.5*10^5;
Y13(i)=ha+tan(pi/180*8)*Xm13(i);
Xm14(i)= i/(numel(z)) * 0.25*10^5;
Y14(i)=ha+tan(pi/180*16)*Xm14(i);
Xm15(i)= i/(numel(z)) * 0.115*10^5;
Y15(i)=ha+tan(pi/180*32)*Xm15(i);
Xm16(i)= i/(numel(z)) * 0.035*10^5;
Y16(i)=ha+tan(pi/180*64)*Xm16(i);
Y17(i)=ha;
end
figure
plot(X,Y,X1,Y1)
figure
plot(X1,Y1,Xm2,Y2,'k',Xm8,Y17,'m')
figure
plot(X1,Y1,'b--',X111,Y111,'b',Xm2,Y2,'k--',Xm22,Y22,'k',Xm8,Y17)     %с учетом ревракции
figure
plot(X111,Y111,'r.-',Xm22,Y22,'k',Xm8,Y17,'m',Xm3,Y3,'g',Xm3,Y33,'b',Xm30,Y30,'g',Xm30,Y330,'b',Xm4,Y4,'g',Xm4,Y44,'b',Xm5,Y5,'g',Xm5,Y55,'b',Xm6,Y6,'g',Xm6,Y66,'b',Xm7,Y7,'g',Xm7,Y77,'b',Xm8,Y8,'g',Xm8,Y88,'b',Xm9,Y9,'g',Xm10,Y10,'g',Xm11,Y11,'g',Xm12,Y12,'g',Xm13,Y13,'g',Xm14,Y14,'g',Xm15,Y15,'g',Xm16,Y16,'g')

%clc;   %очистка командного окна
end