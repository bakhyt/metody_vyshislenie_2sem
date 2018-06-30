function lab0
clc;

a = -1;
b = 0;

eps = 10^(-6);

is_debug = true;
gold_iterations = 1;

fprintf('MINIMUM SEARCH\n');
fplot(@(x) funk(x), [a, b], 'Color', 'b');
hold on;

fprintf('BIT SEARCH METHOD\n');
[x, y] = bitwise_search_minimization(a, b, eps, is_debug, @funk);
p1 = plot(x, y, 'Marker', '.', 'Color', 'y', 'MarkerSize', 15);

fprintf('\nGOLDEN RATIO METHOD\n');
[x, y] = gold_ratio_minimization(a, b, eps, is_debug, @funk);
p2 = plot(x, y, 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);

fprintf('\nPARABOLIC INTERPOLATION METHOD\n');
[x, y] = parabolic_minimization(a, b, eps, is_debug, gold_iterations, @funk);
p3 = plot(x, y, 'Marker', '.', 'Color', 'm', 'MarkerSize', 15);

fprintf('\nNEWTON''S METHOD\n');
[x, y] = newton_minimization(a, b, eps, is_debug, gold_iterations, @funk);
p4 = plot(x, y, 'Marker', '.', 'Color', 'g', 'MarkerSize', 15);

fprintf('\nFMINBIND\n');
[x, fval, ~, output] = fminbnd(@(x) funk(x), a, b);
fprintf('Number of iterations: %d\n', output.iterations);
fprintf('Minimum point: (%7.5f, %7.5f)\n', x, fval);
p5 = plot(x, funk(x), 'Marker', '.', 'Color', 'c', 'MarkerSize', 15);

legend([p1 p2 p3 p4 p5], ...
    {'Bitwise search method', 'Golden section search method', ...
    'Parabolic interpolation method', 'Newton''s method', 'fminbnd'}, ...
    'Location', 'northeast');

hold off;
end

function y = funk(x)
y = sin((-x*x*x*x - 4*x*x*x - 8*x*x -7*x + 1)/sqrt(11)) - 1;
y = y + log10((4*x*x*x*x*x - 4*sqrt(10)*x*x*x*x + 8*x*x*x + 5*x*x - 5*sqrt(10)*x + 9)/(x*x - sqrt(10)*x + 2));
y = -y;
end

% Golden ratio method for initial approximation
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

function f = det(f1, f2, f3, eps)
    f = (-3*f1+4*f2-f3)/2/eps;
end

function f = double_det(f1, f2, f3, eps)
    f = (f3-2*f2+f1)/eps/eps;
end

% Newton's method
function [x, df] = newton_minimization(a, b, eps, is_debug, gold_iterations, func)
[w, e] = golden_iter(a, b, gold_iterations, func);
x = (w + e)/2;

f1 = func(x-eps);
f2 = func(x);
f3 = func(x+eps);
df = det(f1, f2, f3, eps);

iteration = 1;
while (abs(df) > eps)
    if (is_debug)
        fprintf('Iteration no.%d: (%7.5f, %7.5f)\n', iteration, x, f2);
        plot(x, f2, 'Marker', '.', 'Color', 'g', 'MarkerSize', 15);
    end
    iteration = iteration + 1;

    ddf = double_det(f1, f2, f3, eps);
    x = x - df/ddf;
    
    f1 = func(x-eps);
    f2 = func(x);
    f3 = func(x+eps);
    df = det(f1, f2, f3, eps);
end

df = func(x);
fprintf('Minimum point: (%7.5f, %7.5f)\n', x, df);
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
        fprintf('Iteration no.%d: [%7.5f, %7.5f]\n', iteration, x1, x3);
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

fprintf('Minimum point: (%7.5f, %7.5f)\n', x, f);
end

% Golden section search method
function [x, f] = gold_ratio_minimization(a, b, eps, is_debug, func)
t = (sqrt(5) - 1)/2;
l = b - a;

x1 = b - t*l;
x2 = a + t*l;
f1 = func(x1);
f2 = func(x2);

iteration = 1;
while (l > 2*eps)
    if (is_debug)
        fprintf('Iteration no.%d: [%7.5f, %7.5f]\n', iteration, a, b);
    end
    iteration = iteration + 1;
    
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

x = (a + b)/2;
f = func(x);
fprintf('Minimum point: (%7.5f, %7.5f)\n', x, f);
end

% Bitwise search method
function [x, f] = bitwise_search_minimization(a, b, eps, is_debug, func)
s = (b - a)/4;
x0 = a;
f0 = func(x0);

l = true;
iteration = 1;
while (l)
    if (is_debug)
        fprintf('Iteration no.%d: (%7.5f, %7.5f)\n', iteration, x0, f0);
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
fprintf('Minimum point: (%7.5f, %7.5f)\n', x, f);
end