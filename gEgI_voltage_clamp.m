close all
clear
files=dir(fullfile('C:\Users\...','*.abf')); %select folder to analyze pClamp files (voltage clamp data from neurons)

for f=1:length(files)
    eE=-1.5; %calculated reversal potential of AMPA (mV)
    eI=-59.14; %calculated reversal potential of GABA (mV)    
    [channel,smplIntrvl,h]=abfload(files(f).name);
    fullLFP=squeeze(channel(:,1,1:10)); %local field potential (mV)
    fullPatch=squeeze(channel(:,2,1:10)); %intracellular patch clamp (mV)
    fullI=squeeze(channel(:,3,1:10)); %applied current (pA)
   
    
    % uses the transients in holding current to define step on and off
    diffI=diff(fullI); 
    stepOn=find(diffI(:,end)==max(diffI(:,end)));
    stepOff=find(diffI(:,end)==min(diffI(:,end)));
        if stepOn>stepOff
            [stepOff, stepOn] = deal(stepOn, stepOff);
        end
 
    Itheta=round(mean(fullI(stepOn:stepOff,:)));
    Vtheta=[-75:5:(-75+(5*length(Itheta)-1))];
    Patch=fullPatch(stepOn:stepOff,:);
    LFP=fullLFP(stepOn:stepOff,:);
    LFPhilbert=hilbert(LFP);
    LFPphase=angle(LFPhilbert);

    binRange= [-3.2 3.2]; % angle range in radians
 
    edges=linspace(binRange(1),binRange(2),90);
    binAngles=rad2deg(tsmovavg(edges,'s',2));
    binAngles=binAngles(2:end);
    disc=discretize(LFPphase,edges);
    for j=1:length(Patch(1,:));
        for k=1:length(edges)-1;
             phase(k,j)=mean(Patch((disc(:,j)==k),j));
        end
    end

    for i=2:size(phase,1)
        fit=(polyfit(Vtheta,phase(i,:),1))-(polyfit(Vtheta,phase(1,:),1));
        
        eSyn(i)=roots(fit);
        gSyn(i)=fit(1);
        gI(i)=(gSyn(i)*(eE-eSyn(i)))/(eE-eI);
        gE(i)=gSyn(i)-gI(i);
        clear fit
    end

    maxgSyn(f)=max(gSyn);
    gEall(f,:)=gE;
    gIall(f,:)=gI;
    gEPeakAngle(f)=binAngles(find(gE==max(gE)));
    gIPeakAngle(f)=binAngles(gI==max(gI));
    eRev(f)=eSyn(gSyn==max(gSyn));
    gIPeak(f)=gI(gSyn==max(gSyn));
    gEPeak(f)=abs(gE(gSyn==max(gSyn)));
    gRatio(f)=gIPeak(f)./gEPeak(f);
    gRatioAngle(f)=([binAngles(find(gE==max(gE))) binAngles(find(gI==max(gI)))]);

    clearvars -except files gEall gIall binAngles gRatio gRatioAngle gEPeakAngle gIPeakAngle eRev maxgSyn gIPeak gEPeak
    end

for i=1:size(gEall,1)
    for k=1:size(gEall,2)
        gEnorm(i,k)=gEall(i,k)/max(gEall(i,:)+gIall(i,:));
        gInorm(i,k)=gIall(i,k)/max(gEall(i,:)+gIall(i,:));
    end
end


avgExcnorm=median(gEnorm);
avgInhnorm=median(gInorm);
stdExcnorm75=prctile(gEnorm,75);
stdExcnorm25=prctile(gEnorm,25);
stdInhnorm75=prctile(gInorm,75);
stdInhnorm25=prctile(gInorm,25);

avgExc=median(gEall);
avgInh=median(gIall);
meanExc=mean(gEall);
meanInh=mean(gIall);
stdExc=std(gEall);
stdInh=std(gIall);
stdExc75=prctile(gEall,75);
stdExc25=prctile(gEall,25);
stdInh75=prctile(gIall,75);
stdInh25=prctile(gIall,25);


inhIndxnorm=[stdInhnorm75 flip(stdInhnorm25)];
excIndxnorm=[stdExcnorm75 flip(stdExcnorm25)];
binAngles2=[binAngles flip(binAngles)];
figure(1)
% patch(binAngles2,inhIndxnorm,'b','FaceAlpha','.75')

% patch(binAngles2,excIndxnorm,'r','FaceAlpha','.75')
fill(binAngles2,inhIndxnorm,'b')
hold on
fill(binAngles2,excIndxnorm,'r')
plot(binAngles,avgExcnorm,'k','LineWidth',5)
plot(binAngles,avgInhnorm,'k','LineWidth',5)

inhIndx=[stdInh75 flip(stdInh25)];
excIndx=[stdExc75 flip(stdExc25)];
binAngles2=[binAngles flip(binAngles)];
figure(2)
% patch(binAngles2,inhIndx,'b','FaceAlpha','.75')
fill(binAngles2,inhIndx,'b')
hold on
% patch(binAngles2,excIndx,'r','FaceAlpha','.75')
fill(binAngles2,excIndx,'r')
plot(binAngles,avgExc,'k','LineWidth',5)
plot(binAngles,avgInh,'k','LineWidth',5)


