function [x, y] = Circle(r, xc, yc)
    theta = linspace(0,2*pi);
    x = r*cos(theta) + xc;
    y = r*sin(theta) + yc;

end