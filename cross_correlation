clear
close all

files=dir(fullfile('C:\Users\...,'*.abf')); %select folder to analyze abf (pClamp) files
stepSize= 5; %declare size of voltage steps in mV
stepRange= [-110 20]; %declare range limit of voltage steps in mV
voltageStep=(stepRange(1):stepSize:stepRange(2));

for f=1:length(files)

    binsFile=zeros(size(voltageStep)); 

    [channel,smplIntrvl,h]=abfload(files(f).name);

    fullLFP=squeeze(channel(:,1,:));
    fullPatch=squeeze(channel(:,2,:));
    I=squeeze(channel(:,3,:));
    smplFreq=1000000/smplIntrvl;

    dI=diff(I); 
    stepOn=find(dI(:,end)==max(dI(:,end)));
    stepOff=find(dI(:,end)==min(dI(:,end)));
        if stepOn>stepOff
            [stepOff, stepOn] = deal(stepOn, stepOff);
        end

    LFP=fullLFP(stepOn:stepOff,:);
    Patch=fullPatch(stepOn:stepOff,:);
    discCorr=discretize(mean(Patch),voltageStep); % bins voltage response to injected current


    corrTime=0.5*smplFreq;
    
    for i=1:length(LFP(1,:))
        [corr(:,i),lag(:,i)]=xcorr(LFP(:,i)-mean(LFP(:,i)),Patch(:,i)-mean(Patch(:,i)),corrTime,'coeff');
    end
    corrMax=corr(corrTime,:);

    for i=1:length(corrmax)
        binsFile(end+1,discCorr(i))=corrmax(i);
    end

    binsFile(binsFile == 0) = NaN;
    m=nanmean(binsFile);
    maxCorr(f,:)=m;
    clearvars -except maxCorr files edgesCorr
end

hold on

%medAll=nanmedian(binsAll);
%medAll75=prctile(binsAll,75);
%medAll25=prctile(binsAll,25);
%errorbar(edgesCorr,medAll,medAll-medAll25,medAll75-medAll,'capsize',30)
