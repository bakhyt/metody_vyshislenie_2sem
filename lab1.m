function lab1
clc;

a = -1;
b = 0;

eps = 10^(-8);

is_debug = true;

fprintf('MINIMUM SEARCH\n');
fplot(@(x) funk(x), [a, b], 'Color', 'b');
hold on;

fprintf('BITWISE SEARCH METHOD\n');
[x, y] = bitwise_search_minimization(a, b, eps, is_debug, @funk);
p = plot(x, y, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15);
legend((p), {'Bitwise search method'}, 'Location', 'northeast');
   
hold off;
end

% Function and its derivatives (through substraction)

function y = funk(x)
y = sin((x^2 + x - 4)/5);
y = y + cosh((x^3 + 3*(x^2) + 8*x + 8)/(3*x + 9)) - 1.0;
end

function [x, f] = bitwise_search_minimization(a, b, eps, is_debug, func)
s = (b - a)/4;
x0 = a;
f0 = func(x0);

l = true;
iteration = 1;
while (l)
    if (is_debug)
        fprintf('Iteration no.%d: (%7.8f, %7.8f)\n', iteration, x0, f0);
        plot(x0, f0, 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);
    end
    iteration = iteration + 1;
    
    x1 = x0 + s;
    f1 = func(x1);
    if (f0 > f1)
        x0 = x1;
        f0 = f1;
        if (x0 > a && x0 < b)
            continue;
        end
    end
    if (abs(s) <= eps)
        l = false;
    else
        x0 = x1;
        f0 = f1;
        s = -s/4;
    end
end

x = x0;
f = f0;
fprintf('Minimum point: (%7.8f, %7.8f)\n', x, f);
end