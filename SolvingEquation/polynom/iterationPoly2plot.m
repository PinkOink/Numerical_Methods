x=[0.5:0.01:1.0];
y=iterationPoly2(x);
 
figure
plot(x,y)
line([0.0;1.5], [0.0;0.0], 'Color', 'k');
line([0.0;0.0], [0.0;1.0], 'Color', 'k');

line([0.0;1.0], [0.0;1.0], 'Color', 'm');

x1=1.0
line([x1, x1], [0.0, iterationPoly2(x1)],'LineStyle', '--', 'Color', 'r');
line([x1, iterationPoly2(x1)], [iterationPoly2(x1), iterationPoly2(x1)],'LineStyle', '--', 'Color', 'r');
x1 = iterationPoly2(x1)
line([x1, x1], [0.0, iterationPoly2(x1)],'LineStyle', '--', 'Color', 'r');
line([x1, iterationPoly2(x1)], [iterationPoly2(x1), iterationPoly2(x1)],'LineStyle', '--', 'Color', 'r');
text(1.0, 0.98, 'y=x')
text(1.0, 0.6, '\phi(x)')
text(1.0, 0.05, 'x0')
text(x1, 0.05, 'x1')

xlabel('x')
ylabel('y')
grid on
title('Iteration function 2')