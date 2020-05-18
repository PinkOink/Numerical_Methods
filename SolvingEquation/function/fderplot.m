x=[0.1:0.01:50.0];
y=fder(x);

figure
plot(x,y)
line([0;50.0], [0.0;0.0], 'Color', 'k');
line([0.0;0.0], [-1.0;2.0], 'Color', 'k');
xlabel('x')
ylabel('y')
grid on
hold on
title('Derivative')