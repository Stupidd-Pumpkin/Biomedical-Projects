x=0:0.01:10;
y=x.^3;
net=newff(minmax(x),[20,1],{'logsig','purelin','trainlm'});
net.trainparam.epochs=1000;
net.trainparam.goal=1e-25;
net.trainparam.lr=0.01;
net=train(net,x,y);
