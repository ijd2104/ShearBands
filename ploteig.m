%Plor eigen values
hold on
for i = 4:6
    plot(er(i,24:26),ei(i,24:26),'o','MarkerSize',10,'MarkerEdgeColor','g')
    plot(er(i,[24,25,27]),ei(i,[24,25,27]),'*','MarkerSize',10,'MarkerEdgeColor','b')
end
plot([0 0],[-2E10 2E10],'r--')
axis([-2e+03 6e+04 -0.5E10 0.5E10])
legend('Original','Reduced Density')