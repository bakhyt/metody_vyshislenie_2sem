function lab4
clc;

a = -1;
b = 0;

eps = 10^(-2);

is_debug = true;
gold_iterations = 1;

fprintf('MINIMUM SEARCH\n');
fplot(@(x) funk(x), [a, b], 'Color', 'b');
hold on;

fprintf('NEWTON''S METHOD\n');
[x, y] = newton_minimization(a, b, eps, is_debug, gold_iterations, @funk);
p1 = plot(x, y, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15);

fprintf('\nFMINBIND\n');
[x, fval, ~, output] = fminbnd(@(x) funk(x), a, b);
fprintf('Number of iterations: %d\n', output.iterations);
fprintf('Minimum point: (%7.8f, %7.8f)\n', x, fval);
p2 = plot(x, funk(x), 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);
legend([p1 p2], {'Newton''s method', 'fminbnd'}, 'Location', 'northeast');
  
hold off;
end

function y = funk(x)
y = sin((-x*x*x*x - 4*x*x*x - 8*x*x -7*x + 1)/sqrt(11)) - 1;
y = y + log10((4*x*x*x*x*x - 4*sqrt(10)*x*x*x*x + 8*x*x*x + 5*x*x - 5*sqrt(10)*x + 9)/(x*x - sqrt(10)*x + 2));
y = -y;
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

function f = der(f1, f2, f3, eps)
    f = (-3*f1+4*f2-f3)/2/eps;
end

function f = double_der(f1, f2, f3, eps)
    f = (f3-2*f2+f1)/eps/eps;
end

% Newton's method
function [x, df] = newton_minimization(a, b, eps, is_debug, gold_iterations, func)
[w, e] = golden_iter(a, b, gold_iterations, func);
x = (w + e)/2;

f1 = func(x-eps);
f2 = func(x);
f3 = func(x+eps);
df = der(f1, f2, f3, eps);

iteration = 1;
while (abs(df) > eps)
    if (is_debug)
        fprintf('Iteration no.%d: (%7.8f, %7.8f)\n', iteration, x, f2);
        plot(x, f2, 'Marker', '.', 'Color', 'g', 'MarkerSize', 15);
    end
    iteration = iteration + 1;

    ddf = double_der(f1, f2, f3, eps);
    x = x - df/ddf;
    
    f1 = func(x-eps);
    f2 = func(x);
    f3 = func(x+eps);
    df = der(f1, f2, f3, eps);
end

df = func(x);
fprintf('Minimum point: (%7.8f, %7.8f)\n', x, df);
end