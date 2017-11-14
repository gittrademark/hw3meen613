function n = MM(t,y,L1,L2,M,g,w2)
%This function evaluates the value of the Mass Matrix.
n1=[1 0 0 0];
n2=[0 M*L1^2 0 (M*L2*L1/2)*cos(y(1)-y(3))];
n3=[0 0 1 0];
n4=[0 (M*L1*L2)/2*cos(y(3)-y(1)) 0 M*L2^2/3];
%Putting all calculated n's into a single matrix
n=[n1;n2;n3;n4];
end
