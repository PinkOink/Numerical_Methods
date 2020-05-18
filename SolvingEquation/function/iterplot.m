x=[0.0:0.01:5.0];
y=iter(x);

figure
plot(x,y)
line([0;6.0], [0.0;0.0], 'Color', 'k');
line([0.0;0.0], [0.0;2.0], 'Color', 'k');

line([0.0;5.0], [0.0;5.0], 'Color', 'm');

x1=5.0
line([x1, x1], [0.0, iter(x1)],'LineStyle', '--', 'Color', 'r');
line([x1, iter(x1)], [iter(x1), iter(x1)],'LineStyle', '--', 'Color', 'r');
x1 = iter(x1)
line([x1, x1], [x1, iter(x1)],'LineStyle', '--', 'Color', 'r');
line([x1, iter(x1)], [iter(x1), iter(x1)],'LineStyle', '--', 'Color', 'r');

xlabel('x')
ylabel('y')
grid on
hold on
title('Function-iteration')