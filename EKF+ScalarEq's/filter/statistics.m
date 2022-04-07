N = 1000;
str = ['q1','q2','q3','wx','wy','wz'];
j = 1;
for i = 1:6
    varL(i) = var(L(i,:));
    meanL(i) = mean(L(i,:));
    x = min(L(i,:)):(max(L(i,:))-min(L(i,:)))/(N-1):max(L(i,:));
    f = 1/sqrt(2*pi*varL(i))*exp(-(x-meanL(i)).^2/(2*varL(i)));

    [counts1, binCenters1] = hist(L(i,:));
    subplot(2,3,i); bar(binCenters1, counts1,'FaceColor','#4DBEEE');
    hold on;
    plot(x,f*max(counts1)/max(f),'LineWidth', 2, 'Color', '#A2142F'); xlabel(str(j:j+1))
    grid; hold off;
    j = j+2;
end
tak = 1.6602;
delta = sqrt(varL)*tak;
