function lab2
clc;

a = -1;
b = 0;

eps = 10^(-8);

is_debug = true;

fprintf('MINIMUM SEARCH\n');
fplot(@(x) funk(x), [a, b], 'Color', 'b');
hold on;

fprintf('GOLDEN SECTION SEARCH METHOD\n');
[x, y] = gold_ratio_minimization(a, b, eps, is_debug, @funk);
p = plot(x, y, 'Marker', '.', 'Color', 'k', 'MarkerSize', 15);
legend((p), {'Golden section search method'}, 'Location', 'northeast');
   
hold off;
end

function y = funk(x)
y = sin((-x*x*x*x - 4*x*x*x - 8*x*x -7*x + 1)/sqrt(11)) - 1;
y = y + log10((4*x*x*x*x*x - 4*sqrt(10)*x*x*x*x + 8*x*x*x + 5*x*x - 5*sqrt(10)*x + 9)/(x*x - sqrt(10)*x + 2));
y = -y;
end

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
    fprintf('Iteration no.%d: [%7.8f, %7.8f]\n', iteration, x1, x2);
    plot([x1, x2] ,[f1,f2], 'Marker', '.', 'Color', 'r', 'MarkerSize', 15);
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
fprintf('Minimum point: (%7.8f, %7.8f)\n', x, f);
end