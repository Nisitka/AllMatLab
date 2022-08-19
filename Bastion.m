T = 15000;       %весь итервал времени 
to = T/2;      %время в которое произошел отказ
u = 7;    %процент от интервала времени теста, при котором отказы будут синергировать
q = T * u/(100 * 3);

for i=1:T
    f(i) = exp(-(i - to)^2 / (2*q^2) );
end
figure;
plot(1:T,f);

%T = 20000;  %интервал времени теста 
N = 20000;   %кол-во тестов (что эквиваленто запускам или обращениям к функциональному объекту)
n = 50;    %кол-во отказов 
lamda_min = 0.00001;
for i=1:n
   tn(i) = rand*T;  
end

for i=1:T
    m(i) = 0;
    for j=1:n
        m(i) = m(i) + n/N * exp(-(i - tn(j))^2 / (2*q^2));
    end
end
figure
numel(m)
plot(1:T,m,tn,0,'ro');


%вычесление производной
h = 1;
for i=1:(numel(m) - 1*h)
    lamda(i)= (m(i)-m(i+h))/(h);
    if lamda(i) < 0
       lamda(i) = -lamda(i); 
    end
    if lamda(i) < lamda_min
       lamda(i) = lamda_min;
    end
end
figure
plot(1:numel(lamda), lamda);

r = round(numel(lamda)*0.03);

for i=1:(numel(lamda) - 0)
    if (i <= r)
       s = 0;
       for j=1:2*r
          s = s + lamda(i - 1 + j);  
       end
       del = 2*r;
       lamda_(i) = s/del;
    end
    if (i >= numel(lamda) - r)
       s = 0;
       for j=1:2*r
          s = s + lamda(i + 1 - j);  
       end
       del = 2*r;
       lamda_(i) = s/del;
    end
    if (i > r && i < numel(lamda) - r)
        s = 0;
        for j=(i-r):(i+r)
            s = s + lamda(j);
        end
        del = 2*r;
        lamda_(i) = s/del;
    end
end
figure 
plot(1:numel(lamda_), lamda_);

for i=1:numel(lamda_)         
   p(i) = exp(-i*lamda_(i));  
end
figure
plot(1:numel(lamda_), p);


