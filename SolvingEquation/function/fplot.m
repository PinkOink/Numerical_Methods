x=[0.5:0.01:35];
y=f(x);

figure
plot(x,y)
line([0.0;-1.5], [0.0;0.0], 'Color', 'k');
line([-0.5;-0.5], [-1.5;1.5], 'Color', 'k');
xlabel('x')
ylabel('y')
grid on
hold on
title('\phi''(x)')
x0=fzero(@f, 0.5)