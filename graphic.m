function graphic
X = 0:0.001:10;

for i=1:(10/0.001 + 1)
Y(1, i) = exp(X(i));
Y(2, i) = sin(X(i));
end
plot(X,Y(2,:))
end