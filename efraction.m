function efraction
dndz = -6*10^-8;
Rz = 6370000;
Rz2 = Rz/(1 + Rz*dndz);
Ha = 16;
Hb = [50, 100, 500, 1000];
for i=1:4
    Dv1(i)=sqrt(2*Rz2)*(sqrt(Ha)+sqrt(Hb(i)));
    Dv2(i)=sqrt(2*Rz)*(sqrt(Ha)+sqrt(Hb(i)));
end;
plot(Hb,Dv1,Hb,Dv2)
end