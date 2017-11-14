% In order to solve this type of problem, the odefile needs % to
%have switch/case programming with a flag case of %
%'mass'. when the flag is set to default the function % FF1.m
%is called. Later in the time step, the function % with the
%flag of 'mass' is called. As you can see in the % function
%indmot_ode1, when the flag is 'mass', the % function MM1.m
%is called.
function varargout=indmot_ode1(t,y,flag,m,l1,l2)
switch flag
case '' %no input flag
 varargout{1}=FF1(t,y,m,l1,l2);
case 'mass' %flag of mass calls MM1.m
 varargout{1}=MM1(t,y,m,l1,l2);
otherwise
 error(['unknown flag ''' flag '''.']); 
end