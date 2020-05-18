x=[-1.2:0.01:-0.75];
y=iterationPoly1(x);
 
figure
plot(x,y)

line([-0.6;-1.3], [-0.6;-1.3], 'Color', 'm');

x1=-0.75
line([x1, x1], [-0.6, iterationPoly1(x1)],'LineStyle', '--', 'Color', 'r');
text(x1,-0.61,'x0')
line([x1, iterationPoly1(x1)], [iterationPoly1(x1), iterationPoly1(x1)],'LineStyle', '--', 'Color', 'r');
x1 = iterationPoly1(x1)
text(x1,-0.61,'x1')
line([x1, x1], [-0.6, iterationPoly1(x1)],'LineStyle', '--', 'Color', 'r');
line([x1, iterationPoly1(x1)], [iterationPoly1(x1), iterationPoly1(x1)],'LineStyle', '--', 'Color', 'r');
line([iterationPoly1(x1), iterationPoly1(x1)], [-0.6, iterationPoly1(iterationPoly1(x1))],'LineStyle', '--', 'Color', 'r');
line([iterationPoly1(x1), iterationPoly1(x1)], [-0.6, iterationPoly1(iterationPoly1(x1))],'LineStyle', '--', 'Color', 'r');
text(iterationPoly1(x1),-0.61,'x2')
text(-1.2,-1.2,'y=x')
text(-1.24,-0.905,'\phi(x)')
line([-1,-1],[-0.6,-1],'LineStyle', '--', 'Color', 'r');
text(-1,-0.61,'x*')

xlabel('x')
ylabel('y')
grid on
title('Iteration function 1')