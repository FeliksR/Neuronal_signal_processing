close all
clear
files=dir(fullfile('C:\Users\...','*.abf')); %select folder to analyze

for f=1:length(files)

    [channel,smplIntl,h]=abfload(files(f).name);   
    LFP=squeeze(channel(:,1,:)); % set recording channel 1 to local field potential
    Patch=squeeze(channel(:,2,:)); %set recording channel 2 to membrane voltage
    I=squeeze(channel(:,3,:)); % set recording channel 3 to applied current
    samplFreq=1000000/smplIntl; % sampling frequency calculated from sampling interval

    % uses the transients in holding current to define step on and off
    dI=diff(I); 
    stepOn=find(dI(:,end)==max(dI(:,end)));
    stepOff=find(dI(:,end)==min(dI(:,end)));
    if stepOn>stepOff
        [stepOff, stepOn] = deal(stepOn, stepOff);
    end
    
    % Hilbert transform to determine phase of local field potential
    LFPhilbert=hilbert(LFP); 
    LFPphase=angle(LFPhilbert); 
    LFPdiff=diff(LFPphase);
    
    % determines edges of local field potential oscillations in order to bin single cycles
    cycleEdges=[];
    for j=1:length(LFPdiff(1,:))
        for i=1:length(LFPdiff(:,1))
            if LFPdiff(i,j)<-2
                cycleEdges(end+1,j)=i;
            end
        end
    end

    cycle=[];
    for j=1:length(cycleEdges(1,:))
        for i=1:length(cycleEdges(:,1))-1
            if cycleEdges(i+1,j)-cycleEdges(i,j)>(0.1*samplFreq)
                cycle(end+1,j)=cycleEdges(i,j);
                cycle(end+1,j)=cycleEdges(i+1,j);
            end
        end
    end

    for j=1:length(cycle(1,:))
        for i=1:length(cycle(:,1))-1
            if cycle(i,j)>0 && cycle(i+1,j)>0 && mod(i,2)==1
                cycleLength=cycle(i+1,j)-cycle(i,j);
                cycleStart=cycle(i,j);
                cycleEnd=cycle(i+1,j);
                troughHalfWidth=round(cycleLength*0.2); 
                troughStart=cycleStart+troughHalfWidth;
                troughVm(j,cycleStart:troughStart)=Patch((cycleStart:troughStart),j);

            end
        end
    end
    
    troughVm(troughVm==0)=NaN;
    Vm=nanmedian(troughVm);

    stepVm=nanmean(Vm(stepOn+.1*samplFreq:stepOff));
    baseVm=nanmean(Vm(1:stepOn));
    stepI=round(mean(I(stepOn:stepOff,1)));
    imp=1000*((stepVm-baseVm)/(stepI));
    baseVmAll(f)=baseVm;
    stepVmAll(f)=stepVm;
    impAll(f)=imp;
    fexp = @(p,x) p(1) + p(2)*( - exp(-(x)/p(3)));
    p = nlinfit(10:0.1*samplFreq,Vm(stepOn+10:stepOn+0.1*samplFreq),fexp,[1 0 1]);
    tau=p(3)/samplFreq;
    tauAll(f)=tau*1000;

    clearvars -except impAll tauAll files baseVmAll stepVmAll
end

% change_prednqx=[baseVmAll_prednqx; stepVmAll_prednqx]';
% impMean_predqnx=mean(impAll_prednqx);
% impSTD_prednqx=std(impAll_prednqx);
% tauMean_prednqx=mean(tauAll_prednqx);
% tauSTD_prednqx=std(tauAll_prednqx);


% baseVmMean=mean(baseVmAll);
% stepVmMean=mean(stepVmAll);
% baseVmStd=std(baseVmAll);
% stepVmStd=std(stepVmAll);

% boxplot(baseVmAll)
% hold on
% for ii=1:length(baseVmAll)
%   plot([1,2],[baseVmAll(ii),stepVmAll(ii)],'-or',...
%        'MarkerFaceColor',[1,0.5,0.5])
% end
