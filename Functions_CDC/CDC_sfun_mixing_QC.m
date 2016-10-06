function [sys,x0,str,ts] = CDC_sfun_mixing_QC(t,x,u,flag,n_comp,tstep)

% Combination of 2 elements
% In/Output in Q and Concentrations (m3/s, g/m3) as vector  
% sample time: tstep
% vector	Q=1  C= 1..n_comp


switch flag,
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes(n_comp,tstep);
    case 3,
        sys = mdlOutputs(t,x,u,n_comp,tstep);
    case 4,
        sys = -1;
    case {1,2,4,9}
        sys = []; % do nothing
    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end


%=======================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=======================================================================
function [sys,x0,str,ts] = mdlInitializeSizes(n_comp,tstep)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = n_comp+1;
sizes.NumInputs      = 2*(n_comp+1);
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

x0  = [];
str = [];
ts  = [-1 0]; 


%=======================================================================
% mdlOutputs
% Return Return the output vector for the S-function
%=======================================================================
function sys = mdlOutputs(t,x,u,n_comp,tstep)

sys(1) = u(1)+u(n_comp+2);

if sys(1)>0
    sys(2:n_comp+1)=(u(2:n_comp+1)*u(1)+u(n_comp+3:end)*u(n_comp+2))/sys(1);
else

    sys=zeros(n_comp+1,1);
    sys(1) = u(1)+u(n_comp+2);
    if(sys(1)<0)
        warning('negative flows');
        sys(1) = 0;
    end  
end

