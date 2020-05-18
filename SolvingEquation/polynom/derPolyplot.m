x=[-1.25:0.01:1];
y=derPoly(x);
 
figure
plot(x,y)
line([-1.25;1.0], [0.0;0.0], 'Color', 'k');
line([0.0;0.0], [-40.0;50.0], 'Color', 'k');
xlabel('x')
ylabel('y')
grid on
title('Derivative')