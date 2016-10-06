%=============================================================================
% File:			CD1_sfun_Muskingum_Catchment.m
% Purpose:		catchment-transport-function
% Author:		kat
% Date:			050901
% Version		1
%=============================================================================

function [sys,x0,str,ts] = CDC_sfun_Muskingum_Catchment(t,x,u,flag,tstep,CX,CY,n_comp,N)

% Parameters / All provided by block--> IMPORTANT
% tstep     ...Global (discrete) time step for sampling
% K         ...Muskingum Constant [s] (applied to one subreach)
% X         ...Muskingum Constant [-]
% YA,YB     ...Muskingum Constants [ ?? ]
% n=n_comp  ...Number of pollutant components / substances [-]

% N         ...Number of subreahces
N=round(N);

% Internal Variables (Input, State and Output variables)

% x(1)          ...V, Volume in reach [m³]
% x(2..n)       ...C, Component concentration in reach [g/m³]

% u(1)          ...QI, inflow [m³/s]
% u(2..n)       ...CI, Component concentration in inflow [g/m³]

% y(1)          ...QE, outflow [m³/s]
% y(2..n+1)     ...CE, Component concentrations in outflow [g/m³]
% y(n+2)        ...V, Volume in reach [m³]
% y(n+3..2(n+1))...C, Component concentration in reach [g/m³]

switch flag,
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes(n_comp,N,tstep);

    case 2,
        sys = []; % mdlUpdate(t,x,u);

    case 3,
        sys = mdlOutputs(t,x,u,tstep,CX,CY,n_comp,N);

    case 9,
        sys = []; % do nothing

    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end


%=======================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=======================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes(n_comp,N,tstep)

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0; % Nothing, will be integrated outside the block
sizes.NumOutputs     = (1+N)*(n_comp+1); % Volumes and Flows including corresponding concnetrations
%sizes.NumOutputs     = (2)*(n_comp+1);
sizes.NumInputs      = (2+N)*(n_comp+1); % Inflow and volume + concentrations
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = zeros(sizes.NumDiscStates,1);
str = [];
ts  = [-1 0]; % driven by global timesteps defined

% Initial setting of values for previous time steps

u_dat.volume=zeros(N*(n_comp+1));
set_param(gcb,'UserData',u_dat);

% mdlUpdate

function sys = mdlUpdate(t,x,u)

u_dat=get_param(gcb,'UserData');
sys = u_dat.volume;

% mdlOutputs

function sys = mdlOutputs(t,x,u,tstep,CX,CY,n_comp,N)

QI=[];
QL=[];
QE=[];
V =[];
xnew=[];
rates = [];
QL=u(1:(n_comp+1));
QL(1) = QL(1)/N;
QI=u((n_comp+2):(2*n_comp+2));

if(QL(1)<0 || QI(1)<0)
    warning('Negative flows')
end

if(any(QL(2:end)<0) || any(QI(2:end)<0))
    warning('Negative concentrations')
end

x = u(2*(n_comp+1)+1:end);
for j=1:N

    % HYDRAULICS:
    % recalling Vi-1 from previous time step [Vi C1,i C2,i C3,i ... Cn_comp,i]
    Vold=[];
    for m=((j-1)*(1+n_comp)+1):(j*(1+n_comp))
        Vold=[Vold x(m)];
    end
    Vold = Vold';
    % Outflow for current time step
    QE(1)=(QI(1)+QL(1))*CX + Vold(1)*CY;
    % Vi for current time step
%     V(1)=(QI(1)+QL(1)-QE(1))*tstep+Vold(1);
    
    qOut = QE;
    v = Vold(1);
    c = Vold(2:end);
    qIn = QI(1)+ QL(1);
    
    %% Debugging
    if(qIn(1)<0)
        %warning('negative inflow, correcting to zero')
        %qIn(1) = 0;
        
        %  q & c to zero: stable c, noisy q
        %  q to zero: stable c, noisy q
        %  c to zero: unstable c, stable q        
        %  no correction: unstable c, stable q
        %qIn(1) = 0;
%        qIn(2:end) = 0;

%     end
%     if(QI(1)<0)
%         QI(1:end) =0 ;
%     end
%      if(QL(1)<0)
%         QL(1:end) =0 ;
    end   
        %%
        
        
    if(QI(1)+QL(1)~=0)
        cIn = (QI(1)*QI(2:end) + QL(1)*QL(2:end))/(QI(1)+ QL(1));
    else
        cIn = zeros(size(QI(2:end)));
    end
    
    vDot = qIn - qOut;      
        
    %%
    if(v~=0)
        cDot = (qIn*cIn-qOut*c-vDot*c)./v;
    else
        cDot = zeros(n_comp, 1);
    end
    % Calculation of substance flows

    %%
    % Writing new discrete states
    rates =[rates; vDot; cDot];
    
    % Replacing QI with QE for next subreach calculated
    % where outflow of chamber is inflow for the subsequent chamber
  
    QI = [QE;c];
    if(any(isnan(cDot)))
        disp('Nans');
        keyboard
    end
    
    if(c(1)<0.95 && t>4e4)
%           keyboard
    end

end

% renewing the storage vector u_dat for next step usage
% u_dat.volume=rates;
% set_param(gcb,'UserData',u_dat);

% Generating output for block
sys = [QI; rates];



