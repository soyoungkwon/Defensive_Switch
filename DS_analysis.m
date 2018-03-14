
% plot
figure;
bar(dotMatrix(:,1))
hold on;
bar(((dotMatrix(:,3)-dotMatrix(:,5)==0)&(dotMatrix(:,4)-dotMatrix(:,6)==0)).*dotMatrix(:,1),'r')

figure;
bar(dotMatrix(:,1)*6); colormap([.7 .7 .7])
hold on;
bar(((dotMatrix(:,3)-dotMatrix(:,5)==0)+(dotMatrix(:,4)-dotMatrix(:,6))).*dotMatrix(:,1),'b')
bar(((dotMatrix(:,3)-dotMatrix(:,5)==0)&(dotMatrix(:,4)-dotMatrix(:,6)==0)).*dotMatrix(:,1)*6,'r')


figure; bar([mean(dotMatrix(dotMatrix(:,1)==1,end)) mean(dotMatrix(dotMatrix(:,1)==-1,end))])

% escaped/attack
sum(dotMatrix(dotMatrix(:,1)==1,end))
sum(dotMatrix(dotMatrix(:,1)==-1,end))