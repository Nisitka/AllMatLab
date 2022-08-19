function traectoria
digits(100) %точность расчетов 
H0=325; %м  высота полетв ВО



                          %задача 1
%начальные координаты цели (в метрах)
X0=440 * 1000; 
Y0=-237 * 1000;
%перевод в другую систему координат 
k1=[X0; Y0];
M=[cos(pi/2), sin(pi/2); -sin(pi/2), cos(pi/2)];
k2=inv(M)*k1;

V=840 * 1000 / 3600;        %скорость воздушного объекта (м/с)
Qv=1.3;                     %СКО скорости (м/с)
Da=509 * 1000;              %дальность до манерва (в метрах)
T=1;
%k2                         %вывод преобразуемых координат (в метрах)


                                             %задача 2
%формирование точек прямо-го участка до манерва
gh=208 * pi/180;            %курс движения ВО (в радианах)
%--создане массива случайных чисел скорости V c СКО Qv
N=2200;                     %кол-во случайных значений скорости

for i=1:N
    Vn(i)=V + rand*Qv -  rand*Qv;
    if i<=500
    Vn_plot(i) = Vn(i);
    end
end              
figure   
plot(Vn_plot) %вывод графика случайных значений скорости

x(1)=k2(1,1);
y(1)=k2(2,1);
D(1)=0;

i=1;
while D(i)<Da
    x(i+1)=x(i)+Vn(i)*sin(gh)*T;
    y(i+1)=y(i)+Vn(i)*cos(gh)*T;
    D(i+1)=sqrt((x(i+1)-x(1))^2+(y(i+1)-y(1))^2);
    
    i=i+1;
end
%D(i-1)  %вывод конечного значения дальности до манерва
%отображение радиуса работы РЛС
Ox=0;
Oy=0; %координаты РЛС (начало отсчета)
r = 5.5256e+04; %радиус обзора РЛС на данной высоте (нашел в ручную)
rx(1) = -r;

[rx, ry] = Circle(r, 0, 0); %создание массива зоны обнаружения

%x(i)%вывод конечных координат после 1-го прямолинейного участка
%y(i)
figure
plot(x,y,Ox,Oy,'ro',x(1),y(1),'mx',x(i-1),y(i-1),'mx',rx,ry,'r--') %отображение точек прямо-го участка до манерва


                                      %задача 3
%--входные данные 

n=1.42;
left=false; %тип разворота

ghm=169 * pi/180;                  %глубина разворота(в радианах)

Rm= V^2 / (9.8 * sqrt(n^2-1));     %радиус виража (манерва)
%Rm
if left %вычисление центра виража в зависимости от его типа
    xc=x(i)-abs(Rm)*cos(gh);
    yc=y(i)+abs(Rm)*sin(gh);
else
    xc=x(i)-abs(Rm)*cos(gh)*(-1);
    yc=y(i)+abs(Rm)*sin(gh)*(-1);
end
%xc
%yc
figure
plot(x,y,Ox,Oy,'ro',x(1),y(1),'mx',x(i),y(i),'mx',rx,ry,'r--',xc,yc,'bo')

if left %вычис-е нач-го угла движения (в рад) в зав-ти от типа виража
    Fn = 90 * pi/180 + gh;
else
    Fn = 270 * pi/180 + gh; 
end
%Fn/(pi/180) %- 360 %вывод начального значения угла движения в градусах 
Nm = ceil(Rm*ghm/(V*T));       %кол-во точек моделирования разворота
dgh= ghm/Nm;                   %изменения угла для каждого шага
%dgh/(pi/180)           %вывод угла в градусах

%нахождение точек манерва 
F(1)=Fn; 
xm(1)=xc + Rm*sin(F(1));
ym(1)=yc + Rm*cos(F(1));
for j=2:Nm
    if left 
        F(j)=F(j-1)+abs(dgh)*(-1);
    else
        F(j)=F(j-1)+abs(dgh);
    end
    xm(j)=xc + Rm*sin(F(j));
    ym(j)=yc + Rm*cos(F(j));
end
%xm(j)
%ym(j)
figure
plot(x,y,Ox,Oy,'ro',x(1),y(1),'k*',x(i),y(i),'k*',rx,ry,'r--',xc,yc,'ko',xm,ym,'b.-',xm(j),ym(j),'k*')
%начальное значения угла угла 2-го прямоли-го участка
if left 
    gh0=gh-ghm;
else
    gh0=gh+ghm;
end
%gh0/(pi/180)-360 %вывод этого угла в градусах

%получение точек для 2-го прямоли-го участка
x2(1)=xm(j);
y2(1)=ym(j);
D(1)=0;
Db=255*1000;                 %продолжительность 2-го прямолинейного участка
k=1;
while D(k)<Db
    x2(k+1)=x2(k)+Vn(k)*sin(gh0)*T;
    y2(k+1)=y2(k)+Vn(k)*cos(gh0)*T;
    D(k+1)=sqrt((x2(k+1)-x2(1))^2+(y2(k+1)-y2(1))^2);

    k=k+1;
end
k=k-1;
%x2(k)
%y2(k)
%D(k-1) %вывод конечного значения дальности 2-го участка
figure
plot(x,y,'b',Ox,Oy,'ro',x(1),y(1),'k*',xm,ym,'b',x2,y2,'b',x2(k),y2(k),'k*')

%объединение всех точек траектории палета ВО
for p=1:i
    X(p) = x(p);
    Y(p) = y(p);
end
for p=i+1:(i+j)
    X(p) = xm(p-i);
    Y(p) = ym(p-i);
end
for p=i+j+1:(i+j+k)
    X(p) = x2(p-(i+j));
    Y(p) = y2(p-(i+j));
end
%plot(X,Y)    %проверка траектории (вся идеальная траектория)
                       %10. Моделирование процесса обнаружения (глава 2)
                     
R = 55256;                 %знач. абциссы в точке пересечения траектории полета ВО и ЗО (в метрах)                       
dzo = sqrt(R^2 + H0^2)      %дальность до этой точки, т.е. дальность ЗО
Nshagov = p;                 %сколько вообще шагов в идеальной траектории

k = 1;
N_A0 = 0;
N_A1 = 0;


for i=1:10:Nshagov  %каждую десятую смоделированную точку
    r(k)=sqrt(X(i)^2 + Y(i)^2 + H0^2);  %дальность до ВО
    if r(k)<=dzo            %входит ли ВО в ЗО? 
        if k==1 
            ko = i;        % запоминаем на каком шаге ВО вошел в ЗО впервые
        end
        P(k)=exp(-0.68*(r(k)/dzo)^4);   %вероятность обнаружения ВО
        Tr(k, 1) = X(i);
        Tr(k, 2) = Y(i);
        Tr(k, 3) = P(k);
        if P(k)>= rand                  %модулирование процесса обнаружения
            Tr(k, 4) = 1;
            N_A1 = N_A1 + 1;     %подсчет кол-ва обнаружений
        else
            Tr(k, 4) = 0;
            N_A0 = N_A0 + 1;
        end
       
        if r(k)>dzo
           
            break;
        end
        k=k+1;
    end
end  
Nzo = k - 1;  %кол-во точек в ЗО (совпадает с N_A0 + N_A1)
clear k;
%N_A1       - кол-во обнаружений
%N_A0       - кол-во не обнаружений
%figure
%plot(1:Nzo,P)

%выявление точек обнаружения и запись их:
s1=0;
s2=0;
for i=1:Nzo
    if Tr(i,4) == 1
        s1 = s1 + 1;
        X_A1(s1) = Tr(i,1);
        Y_A1(s1) = Tr(i,2);
        
    else
        s2 = s2 + 1;
        X_A0(s2) = Tr(i,1);
        Y_A0(s2) = Tr(i,2);
    end
    
end
clear s1;
clear s2;
figure
plot(X,Y,'r--',X_A1,Y_A1,'go',X_A0,Y_A0,'ro',Ox,Oy,'ko',rx,ry,'b--')
%                              11. Моделирование процесса измерения в РЛС
%                              координат воздушного объекта в сверической
%                              системе координат

c = 3 * 10^8; %м/с
Tu = 5 / 10^6; %с - время импульса
Qp = 4.6; %ширина диаграммы направленности (град) 
Emin = 0.1;    %минимальный угол места
%ошибки измерений:
Qr = (c * Tu) / (4 * sqrt(3)) %СКО наклонной дальности
Qaz = 0.2 * Qp                %СКО азимута (град)
Qym = 0.1 * Qp                %СКО угла места (град)
%забиваю уже умеющиеся значения в матрицу Izm:
%Nzo
%numel(r)
%ko
for i=1:Nzo
    %записываем координыты ВО в ЗО
    Izm(i, 1) = X(ko + (i-1)*10);  
    Izm(i, 2) = Y(ko + (i-1)*10); 
    
    Izm(i, 3) = r(i);                         %r истинное 
    Izm(i, 4) = r(i) + rand * Qr - rand * Qr; %r с крышкой
    rrr(i) = Izm(i, 4) - r(i);   %разница r и r^

    Izm(i, 5) = asin(H0/Izm(i, 3)) * 180/pi;
    Izm(i, 6) = Izm(i, 5) + rand * Qym - rand * Qym; %E^
    if Izm(i, 6)<= 0
        Izm(i, 6)= Emin;
    end
    %истенное значение азимута:
    if  Izm(i, 1) > 0 & Izm(i, 2) > 0
        Izm(i, 7) = atan(Izm(i, 1)/Izm(i, 2)) *180/pi;
    end
    if  Izm(i, 1) > 0 & Izm(i, 2) < 0
        Izm(i, 7) = 90 + atan(abs(Izm(i, 2))/Izm(i, 1)) *180/pi;
    end
    if  Izm(i, 1) < 0 & Izm(i, 2) < 0
        Izm(i, 7) = 180 + atan(Izm(i, 1)/Izm(i, 2)) *180/pi;
    end
    if  Izm(i, 1) < 0 & Izm(i, 2) > 0
        Izm(i, 7) = 270 + atan(Izm(i, 2)/abs(Izm(i, 1))) *180/pi;
    end
    Izm(i, 8) = Izm(i, 7) + rand*Qaz - rand*Qaz; %азимут с крышкой
    Izm(i, 9) = Tr(i, 4);       %A

end
figure
plot(Izm(:,1),Izm(:,2),'ro')
figure
plot(1:Nzo,Izm(:, 3),'r',1:Nzo,Izm(:, 4),'bo');       %истинная дальность до цели
figure
plot(1:Nzo,rrr,'r'); 
figure
plot(1:Nzo,Izm(:,5),'r',1:Nzo,Izm(:,6),'bo')    %угол места
figure
plot(1:Nzo, Izm(:, 6) - Izm(:, 5),'r');
figure
plot(1:Nzo,Izm(:, 7),'r',1:Nzo,Izm(:, 8),'bo')    %азимут
figure
plot(1:Nzo, Izm(:, 8) - Izm(:, 7),'r');
%plot(1:Nzo, Izm(:, 9));
%                Моделтрование процесса обнаружения траектории
k = 3;    %кол-во обнаружений
m = 4;    %из скольки измерений 
l = 1;
Ntra = 0;

while l<=(Nzo - m + 1) %проход всех точек в ЗО с учетом минимального кол-ва точек для обнаружения
    if Izm(l, 9) == 1                 %если ВО обнаружин, то проверяем на m                                
        S = 0;  %переменная для подсчета кол-ва обнаружений подряд
        for c=1:m
           S = S + Izm(l + c - 1, 9);
        end
        if S >= k                     %если удолитворяет критерию k/m:
            pop = l;                  %запоминаем на каком шаге выполняется условие завязки
            for i=1:(Nzo - l + 1)
                Trob(i, 1) = Izm(l + i - 1, 1); %X
                Trob(i, 2) = Izm(l + i - 1, 2); %Y
                Trob(i, 3) = Izm(l + i - 1, 4); %r^
                Trob(i, 4) = Izm(l + i - 1, 6); %ym^
                Trob(i, 5) = Izm(l + i - 1, 8); %az^
                Trob(i, 6) = Izm(l + i - 1, 9); %A^
            end
            Ntra = i - 1; %Nzo - l + 1;  %сколько точек в завязанной траектории 
            l = Nzo; %чтобы выйти из while
        else
            l = l + 1;
        end
    else
        l = l + 1;
    end
end

clear h;
clear l;
figure
plot(X,Y,'r--',X_A1,Y_A1,'g*',Ox,Oy,'ro',rx,ry,'b--',X_A0,Y_A0,'r*',Trob(:, 1),Trob(:, 2),'bs')

%             Пересчет координат воздушного объекта из системы координат
%             РЛС в систему координат автоматизированного пункта управления
k = 0;
for j=1:Ntra
   if Trob(j, 6) == 1 %ВО обнаружен? 
        k = k + 1;
        r_gor(k) =  Trob(j, 3) * cos(Trob(j, 4)* pi/180);
        r_gor_is(k) = Izm(j+pop-1, 3) * cos(Izm(j+pop-1, 5)* pi/180);
        X_is(k) = r_gor_is(k) * sin(Izm(j+pop-1, 7)* pi/180);
        Y_is(k) = r_gor_is(k) * cos(Izm(j+pop-1, 7)* pi/180);
        Z_is(k) = H0;
        H_is(k) =  Z_is(k) + (Izm(j+pop-1, 3))^2 / (2*12989*10^3);
        
        X_(k) = r_gor(k) * sin(Trob(j, 5)* pi/180);
        Y_(k) = r_gor(k) * cos(Trob(j, 5)* pi/180);
        Z_(k) = Trob(j, 3) * sin(Trob(j, 4)* pi/180);
        H_(k) = Z_(k) + (Trob(j, 3))^2 / (2*12989*10^3);
        if H_(k) < 0
             H_(k) = 0;
        end 
   else
        k = k + 1; 
        
        r_gor(k) =  Trob(j, 3) * cos(Trob(j, 4)* pi/180);
        r_gor_is(k) = Izm(j+pop-1, 3) * cos(Izm(j+pop-1, 5)* pi/180);
        X_is(k) = r_gor_is(k) * sin(Izm(j+pop-1, 7)* pi/180);
        Y_is(k) = r_gor_is(k) * cos(Izm(j+pop-1, 7)* pi/180);
        Z_is(k) = H0;
        H_is(k) =  Z_is(k) + (Izm(j+pop-1, 3))^2 / (2*12989*10^3);
        
        X_(k) = 0;
        Y_(k) = 0;
        Z_(k) = 0;
        H_(k) = 0;                
   end

end
figure
plot(Izm(:,1),Izm(:,2),'r.-',Ox,Oy,'k^',rx,ry,'b--',X_,Y_,'bo')
figure
plot(1:k, X_, 'bo',1:k,X_is,'r.-')
figure
plot(1:k, Y_, 'bo',1:k,Y_is,'r.-')
figure
plot(1:k, Z_, 'bo',1:k,Z_is,'r--',1:k, H_,'k.--')

%                              Преобразование из МПСК РЛС в МПСК ПОИ
X_poi = -30000;
Y_poi = 11000;
Jh = 0;

K0 = [-X_poi; -Y_poi];
M = [cos(Jh), -sin(Jh); sin(Jh), cos(Jh)];
for i=1:k
   
   Kpoi = K0 + M * [X_(i); Y_(i)]; 
   X1_(i) = Kpoi(1,1);
   Y1_(i) = Kpoi(2,1);
   
   Kpoi = K0 + M * [X_is(i); Y_is(i)]; 
   X1_is(i) = Kpoi(1,1);
   Y1_is(i) = Kpoi(2,1);
end
%перевод в другие координаты радиус ЗО
rx = rx - X_poi;
ry = ry - Y_poi;

figure
plot(X1_,Y1_,'.',-X_poi,-Y_poi,'o',rx,ry,'b--',X_poi,Y_poi,'-s','MarkerSize',10,...
   'MarkerFaceColor',[1 .6 .6]);
figure
plot(1:numel(X1_),X1_,'b.:',1:numel(X_),X_,'r--') 
figure
plot(1:numel(Y1_),Y1_,'b.:',1:numel(Y_),Y_,'r--') 

%         14.Экстраполяция пар-ов траектории по фиксированной выборке 
%            измеренных координат
l = 3;   %размер выборки 
for i=1:numel(X1_)
    X1_ex(i) = 0;
    Y1_ex(i) = 0;
    H_ex(i) = 0;
    if i>l 
        for j=1:l
          %координата X
          if X1_(i-(l-j+1)) == (0 - X_poi)
          X1_ex(i) = X1_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * X1_ex(i-(l-j+1)); %если "не обнаружение" то берем экс-ю точку
          else
          X1_ex(i) = X1_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * X1_(i-(l-j+1));    
          end
          %координата Y
          if Y1_(i-(l-j+1)) == (0 - Y_poi)
          Y1_ex(i) = Y1_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * Y1_ex(i-(l-j+1));
          else
          Y1_ex(i) = Y1_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * Y1_(i-(l-j+1));    
          end
          %координата Z
          if H_(i-(l-j+1)) == 0
          H_ex(i) = H_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * H_ex(i-(l-j+1));
          else
          H_ex(i) = H_ex(i) + (6*j - 2*l - 4)/(l*(l-1)) * H_(i-(l-j+1));    
          end
    
        end
    else  %до выборки записывать измереенные координаты
        if X1_(i) == (0 - X_poi)
            X1_ex(i) = X1_(i-1);
        else
            X1_ex(i) = X1_(i);
        end
        
        if Y1_(i) == (0 - Y_poi)
            Y1_ex(i) = Y1_(i-1);
        else
            Y1_ex(i) = Y1_(i);
        end
        
        if X1_(i) == 0
            H_ex(i) = H_(i-1);
        else
            H_ex(i) = H_(i);
        end
    end  
end
figure %координата X
plot(1:numel(X1_),X1_,'bo-',1:numel(X1_ex),X1_ex,'k.-',1:numel(X1_is),X1_is,'r.--')
figure %координата Y
plot(1:numel(Y1_),Y1_,'bo-',1:numel(Y1_ex),Y1_ex,'k.-',1:numel(Y1_is),Y1_is,'r.--')
figure %координата Z
plot(1:numel(H_),H_,'bo-',1:numel(H_ex),H_ex,'k.-',1:numel(H_is),H_is,'r.--')

%замена необнаруженных точек на экстрополированные
for i=1:numel(X1_)
   
   if X1_(i) == (0 - X_poi)
      X1_(i) = X1_ex(i); 
   end
   if Y1_(i) == (0 - Y_poi)
      Y1_(i) = Y1_ex(i); 
   end
   if H_(i) == 0 & H_ex(i)>0
      H_(i) = H_ex(i); 
   end
   
end
figure
plot(1:numel(X1_),X1_,'bo-',1:numel(X1_is),X1_is,'r.--')
figure
plot(1:numel(Y1_),Y1_,'bo-',1:numel(Y1_is),Y1_is,'r.--')
figure
plot(1:numel(H_),H_,'bo-',1:numel(H_is),H_is,'r.--')
figure
plot(X1_,Y1_,'o:k',X1_is,Y1_is,'.',-X_poi,-Y_poi,'ro',rx,ry,'b--',X_poi,Y_poi,'r-s')

%            15.Оценка параметров траекторий по фиксированной выборке
%            измеренных координат 
l = 6;
for i=1:numel(X1_)
    X1_sg(i) = 0;
    Y1_sg(i) = 0;
    H_sg(i) = 0;
    if i>l 
        for j=1:l
  
          X1_sg(i) = X1_sg(i) + (6*j - 2*l - 2)/(l*(l+1)) * X1_(i-(l-j));
          Y1_sg(i) = Y1_sg(i) + (6*j - 2*l - 2)/(l*(l+1)) * Y1_(i-(l-j));
          H_sg(i) = H_sg(i) + (6*j - 2*l - 2)/(l*(l+1)) * H_(i-(l-j));
             
        end
    else   %если размера окна не хватает то принимать незглаженное значение
    X1_sg(i) = X1_(i);
    Y1_sg(i) = Y1_(i);
    H_sg(i) = H_(i);    
    end  
end
figure
plot(1:numel(X1_),X1_,'ko:',1:numel(X1_is),X1_is,'r.-',1:numel(X1_),X1_sg,'sm-.')
figure
plot(1:numel(Y1_),Y1_,'ko:',1:numel(Y1_is),Y1_is,'r.-',1:numel(Y1_),Y1_sg,'sm-.')
figure
plot(1:numel(H_),H_,'ko:',1:numel(H_is),H_is,'r.-',1:numel(H_),H_sg,'sm-.')
figure
plot(X1_,Y1_,'ks:',X1_is,Y1_is,'ro-',-X_poi,-Y_poi,'ro',rx,ry,'b--',X_poi,Y_poi,'r-s',X1_sg,Y1_sg,'m^--')
%                  16.Оценка показателей качества фильтрации дискретных
%                  фазовых координат траекторий ВО
k = 0;
for i=1:numel(Trob(:, 6))
    if Trob(i, 6) == 1   %ждем пока точка окажется в ЗО
        k = k + 1;
        Qx(k) = sqrt((Trob(i, 3) *  Qaz * pi/180)^2);
        Qy(k) = Qr;
        Qh(k) = sqrt((Trob(i, 3) *  Qym * pi/180)^2);
        
        Qx_(k) = sqrt((X1_is(k) - X1_sg(k))^2);
        Qy_(k) = sqrt((Y1_is(k) - Y1_sg(k))^2);
        Qh_(k) = sqrt((H_is(k) - H_sg(k))^2);
        
    end
    
end
Qx_(1)
figure
plot(1:k,Qx,'.b:',1:k,Qx_,'sr:')
figure
plot(1:k,Qy,'.b:',1:k,Qy_,'sr:')
figure
plot(1:k,Qh,'.b:',1:k,Qh_,'sr:')
%                     17. Оценка значения курса полета воздушного объекта
clear gh;
for i=1:(numel(X1_) - 1)
    if X1_(i+1) >= X1_(i)
        gh_(i) = acos(( Y1_(i) - Y1_(i+1)) / sqrt( (X1_(i) - X1_(i+1))^2  +  (Y1_(i) - Y1_(i+1))^2 )) * 180/pi;                                   
        
        
    else
         gh_(i) = 360 - acos(( Y1_(i) - Y1_(i+1)) / sqrt( (X1_(i) - X1_(i+1))^2  +  (Y1_(i) - Y1_(i+1))^2 )) * 180/pi;        
         
    end
    
    if X1_is(i+1) >= X1_is(i)
         gh_is(i) = acos(( Y1_is(i) - Y1_is(i+1)) / sqrt( (X1_is(i) - X1_is(i+1))^2  +  (Y1_is(i) - Y1_is(i+1))^2 )) * 180/pi;
    else
         gh_is(i) = 360 - acos(( Y1_is(i) - Y1_is(i+1)) / sqrt( (X1_is(i) - X1_is(i+1))^2  +  (Y1_is(i) - Y1_is(i+1))^2 )) * 180/pi; 
    end
end
figure
plot(1:numel(gh_), gh_, 'bo:',1:numel(gh_), gh_is, 'r--')

for i=1:numel(gh_)
    d_gh(i) = gh_is(i) - gh_(i);
    
end
figure
plot(1:numel(gh_),d_gh,'ro--');
mat = mean(d_gh)
sko = std(d_gh)
mediana = median(d_gh)
end