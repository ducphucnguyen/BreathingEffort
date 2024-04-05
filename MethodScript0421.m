% Code to generate an attempted flow signal from an oesophageal pressure
% signal
% Requires Flow, Oesophageal pressure (Poes) and Time signals at same sample rate
[Vdot_intended, V_intended, E, R, Residual, BreathTable] = generateAttemptedFlow(Flow, Poes, Time)
%% Fit E and R from stable breathing period
%t1= start index of stable period 
%t2= end index of stable period 
%allow you to manually select stable flow data to fit E and R

figure('Name','Select stable flow and pressure region to fit E and R');
subplot(2,1,1);plot(Flow)
subplot(2,1,2);plot(Poes)
[t,y]=ginput(2);
figure('Name','Refine selected region to fit E and R');
t2=round(t(2))
t1=round(t(1))
subplot(2,1,1);plot(Flow(t1:t2))
subplot(2,1,2);plot(Poes(t1:t2))
[t,y]=ginput(2);
t1=round(t(1))+t1;
t2=round(t(2)+t2);

dTime=Time(2)-Time(1);
Volume=0;
 for i=1:1:round(t2-t2)
   Volume(i+1)= Volume(i) + Flow(i+t1)*dTime;
  end

Volume=detrend(Volume)';
Pmus=Poes(t1:t2)-10*Volume;

tbl=table(Time(t1:t2), Pmus, Flow(t1:t2), Volume,...
     'VariableNames',{'Time','Poes','Flow','Vol'});
 lm=fitlm(tbl,'Poes~Flow+Vol');
 R=(lm.Coefficients.Estimate(2))
 E=(lm.Coefficients.Estimate(3))
 Residual=(lm.Rsquared.Adjusted)
 
%% Generate attempted flow and volume signals over all data using E and R in model

t1=1;
t2=length(Time);

%Apply volume correction term to generate Pmus singal for whole night 
%Using established chest wall compliance correction 10 cmH2O/l 
Volume=0;
 for i=t1:t2-1
   Volume(i+1)= Volume(i) + Flow(i)*dTime;
  end
Volume=detrend(Volume);
Volume=Volume';
y=bandpass(Volume,[0.2 5],1/dTime);
Pmus=Poes-10*y;

%initialize attempted volume and flow signals
V_intended=zeros(length(Time),1);
Vdot_intended=zeros(length(Time),1);
%Loop to simulate intended flow / intended volume traces 
startidx=[];endidx=[]
for i=1:(t2-t1)
    Vdot_intended(i) = (-(-Pmus(i))-(E*V_intended(i)))/R; %intended flow
    if i<length(Flow)
        V_intended(i+1) = V_intended(i) + Vdot_intended(i)*dTime; %intended vol
    end
     if (Vdot_intended(i)/Vdot_intended(i-1))<0 
        if Vdot_intended(i)>0
            startidx=[startidx; i];
        elseif Vdot_intended(i)<0
            endidx=[endidx;i];
        end
    end
end

%Refine attempted effort start and end points
if endidx(1)<startidx(1)
    endidx(1)=[]
end
if length(endidx)<length(startidx)
    startidx(end)=[]
elseif  length(endidx)>length(startidx)
    endidx(end)=[]
end

%% Calculate effort by effort attempted and achieved ventilation, and get flow:effort ratio

BreathTable=table;
BreathTable.InspStartIdx=startidx;
BreathTable.InspEndIdx=endidx;
BreathTable.BreathDuration=(endidx-startidx).*dTime;
BreathTable(BreathTable.BreathDuration<1,:)=[]; %Effort timing threshold: set at minimum duration 1 second
BreathTable(end,:)=[] ;


for jj=1:height(BreathTable)
    try
    data=Flow(BreathTable.InspStartIdx(jj)-20:BreathTable.InspEndIdx(jj)+20);
    BreathTable.VI(jj)=nansum(data(data>0))*dTime/60; %Achieved minute ventilation
    BreathTable.EffortVI(jj)=nansum(Vdot_intended(Breaths(jj,1):Breaths(jj,2)))*dTime/60; %attempted minute ventilation
    catch
    end
end
 
BreathTable.FlowEffortRatio=BreathTable.VI./BreathTable.EffortVI
BreathTable.FlowEffortRatio(BreathTable.FlowEffortRatio>100)=100;
BreathTable.FlowEffortRatio(BreathTable.FlowEffortRatio<0)=0;

 