function yp=FF(t,y,L1,L2,M2,g,w2)
yp=zeros(4,1);
yp(1) = y(2);
yp(2) = -w2*L1*sin(y(1)) - (M2*L1*L2/2)*(y(4)^2)*sin(y(1)-y(3));
yp(3) = y(4);
yp(4) = (-w2*L2*sin(y(3)))/2 + (M2*L1*L2/2)*(y(2)^2)*sin(y(1)-y(3));
end