
N = 6;
for i=1:N
    Ax(i) = 5.0 + 5.0*rand;
    Ay(i) = 5.0 + 5.0*rand;
    Bx(i) = 5.0*rand;
    By(i) = 5.0*rand;
end
figure 
plot(Ax,Ay,'b+',Bx,By,'ro');

Ao_x = 0;
Ao_y = 0;
Bo_x = 0;
Bo_y = 0;
for i=1:N
    Ao_x = Ao_x + Ax(i)/N;
    Ao_y = Ao_y + Ay(i)/N;
    Bo_x = Bo_x + Bx(i)/N;
    Bo_y = Bo_y + By(i)/N;
end;
Line(1,:) =  [Ao_x, Bo_x];
Line(2,:) =  [Ao_y, Bo_y];
figure 
plot(Ax,Ay,'b+',Bx,By,'ro',Line(1,:),Line(2,:),'*k-');  

A = [Ao_x, Ao_y];
B = [Bo_x, Bo_y];
W = A - B;
modA = sqrt(A(1)^2 + A(2)^2);
modB = sqrt(B(1)^2 + B(2)^2);
%t = (modA^2 - modB^2)/2;

k = (A(2)-B(2)) / (A(1)-B(1));
x = 0:0.1:10;
y = W(2) - (x - W(1))/k;

figure
plot(Ax,Ay,'b+',Bx,By,'ro',Line(1,:),Line(2,:),'*k-',x,y,'--'); 

