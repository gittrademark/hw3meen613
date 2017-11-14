function varargout=indmot_ode(t,y,flag,L1,L2,M2,g,w2)

switch flag
case '' %no input flag
varargout{1}=FF(t,y,L1,L2,M2,g,w2);
case 'mass' %flag of mass calls mass.m
varargout{1}=MM(t,y,L1,L2,M2,g,w2);
otherwise
error(['unknown flag ''' flag '''.']);
end