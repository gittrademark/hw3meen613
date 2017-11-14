% the following function contains the mass matrix.
%it is separately stored in a file named, MM1.m
function n = MM1(t,y,m,L1,L2)
n1=[1 0 0 0 0];
n2=[0 m*L1^2 0  (m*L2*L1/2).*cos(y(1)-y(3)) L1.*sin(y(1))];
n3=[0 0 1 0 0];
n4=[0 (m*L1*L2/2).*cos(y(1)-y(3)) 0 (m*L2^2)/3 (L2)*sin(y(3))];
n5=[0 L1*sin(y(1)) 0 L2*sin(y(3)) 0];
n=[n1;n2;n3;n4;n5]; 
end