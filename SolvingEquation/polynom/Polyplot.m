x=[-1.25:0.01:1];
y=Poly(x);
 p=[9.4481, 5.0141, -4.367, -0.2807, -0.3477]
figure
plot(x,y)
line([-1.25;1.0], [0.0;0.0], 'Color', 'k');
line([0.0;0.0], [-2.0;10.0], 'Color', 'k');
xlabel('x')
ylabel('y')
grid on
title('Polynom')
x0=roots(p)