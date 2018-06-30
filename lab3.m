function lab3
clc;

a = -1;
b = 0;

eps = 10^(-6);

is_debug = true;
gold_iterations = 1;

fprintf('MINIMUM SEARCH\n');
fplot(@(x) funk(x), [a, b], 'Color', 'b');
hold on;

fprintf('PARABOLIC INTERPOLATION METHOD\n');
[x, y] = parabolic_minimization(a, b, eps, is_debug, gold_iterations, @funk);
p = plot(x, y, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15);
legend((p), {'Parabolic interpolation method'}, 'Location', 'northeast');
   
hold off;
end

function y = funk(x)
y = sin((x^2 + x - 4)/5);
y = y + cosh((x^3 + 3*(x^2) + 8*x + 8)/(3*x + 9)) - 1.0;
end

function [w, e] = golden_iter(a, b, iterations, func)
t = (sqrt(5) - 1)/2;
l = b - a;

x1 = b - t*l;
x2 = a + t*l;
f1 = func(x1);
f2 = func(x2);

for j = 1:iterations
    if (f1 > f2)
        a = x1;
        l = b - a;
        x1 = x2;
        f1 = f2;
        x2 = a + t*l;
        f2 = func(x2);
    else
        b = x2;
        l = b - a;
        x2 = x1;
        f2 = f1;
        x1 = b - t*l;
        f1 = func(x1);
    end
end

w = a;
e = b;
end

function [x, f] = parabolic_minimization(a, b, eps, is_debug, gold, func)
[w, e] = golden_iter(a, b, gold, func);

x1 = w;
x2 = (w + e)/2;
x3 = e;

f1 = func(x1);
f2 = func(x2);
f3 = func(x3);

a1 = (f2 - f1)/(x2 - x1);
a2 = ((f3 - f1)/(x3 - x1) - (f2 - f1)/(x2 - x1))/(x3 - x2);
x = (x1 + x2 - a1/a2)/2;
f = func(x);

l = true;
iteration = 1;
while (l)
    if (is_debug)
        fprintf('Iteration no.%d: [%7.8f, %7.8f]\n', iteration, x1, x3);
        plot(x2, f2, 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);
    end
    iteration = iteration + 1;
    
    xs = x;
    if (x < x2)
        if (f > f2)
            x1 = x;
            f1 = f;
        else
            x3 = x2;
            f3 = f2;
            x2 = x;
            f2 = f;
        end
    else
        if (f >= f2)
            x3 = x;
            f3 = f;
        else
            x1 = x2;
            f1 = f2;
            x2 = x;
            f2 = f;
        end
    end
    
    a1 = (f2 - f1)/(x2 - x1);
    a2 = ((f3 - f1)/(x3 - x1) - (f2 - f1)/(x2 - x1))/(x3 - x2);
    x = (x1 + x2 - a1/a2)/2;
    f = func(x);
    l = (abs(x - xs) > eps);
end

fprintf('Minimum point: (%7.8f, %7.8f)\n', x, f);
end