%the following function contains the right hand side of the
%differential equation of the form
%M(t,y)*y'=F(t,y)
%i.e. it contains F(t,y).it is also stored in a separate %file
%named, FF1.m.
function yp=FF1(t,y,m,L1,L2)
g=9.81;
w=m*g;
yp=zeros(5,1);
yp(1)=y(2);
yp(2)=-m*L1*L2/2.*sin(y(1)-y(3)).*y(4)^2 - w*L1.*sin(y(1));
yp(3)=y(4);
yp(4)=m*L1*L2/2.*sin(y(1)-y(3)).*y(2)^2 - w*L2/2.*sin(y(3)); %theta
yp(5)= -L1*cos(y(1)).*(y(2).^2)-L2*cos(y(3)).*(y(4).^2); 
end