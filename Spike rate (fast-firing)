close all
clear
files=dir(fullfile('C:\Users\...','*.abf')); %select folder to analyze pClamp file (.abf)

spkISI=[];

for f=1:length(files)
    
    %load data into variables
    [channel,smplIntvl,h] = abfload(files(f).name);
    
    smplFreq = (1000000/smplIntvl);
    
    vExt=squeeze(channel(:,1,:)); % Extracellular spike recording
    temp=squeeze(channel(:,2,:)); % Bath temperature probe
     
    
    %bandpass-filter data 
    filterOrder = 500; 
    bandpassLow = 100; %in Hz
    bandpassHigh = 5000; %in Hz

    filter = fir1(filterOrder, [bandpassLow/(smplFreq/2) bandpassHigh/(smplFreq/2)]);
    vExtFilt = filtfilt(filter, 1, vExt);
    vExtChunk = round(length(vExtFilt)/20);

    %sets spike detection threshold relative to amplitude of signal 
    for i=1:18
        vMax(i) = max(vExtFilt(vExtChunk*(i):vExtChunk*(i+1)));
    end
    thresh=(2/5) * (median(vMax));
    
    % index of spike times
    spkIndx = 0;
    for i=2:length(vExtFilt)-1000
        if i>spkIndx(end) + (.001*smplFreq) && (sum(vExtFilt(i:(i+1)) > thresh)>1)
            spkIndx(end+1) = i;
        end
    end
    spkISI = [spkISI (diff(spkIndx)./smplFreq)];
    spkDiff = (1./(diff(spkIndx)./smplFreq));

    rate(f) = mean(spkDiff); %spike rate
    coefVar(f) = std(diff(spkIndx(2:end)))/mean(diff(spkIndx(2:end))); % spike CV
        
    [ac,xbin] = acf(spkIndx./smplFreq,.001,1);
    acStr(f) = max(ac(1010:length(ac))); %autocorrelation strength

    % to detect signal loss in recordings, plots the anomolous recording for review
    if(coefVar(f))<0.09
        fh = figure(f);
        fh.WindowState = 'maximized';
        plot(diff(spkIndx)./smplFreq);
        fh = figure(f*100);
        fh.WindowState = 'maximized';
        plot(linspace(1,length(vExt)/smplFreq,length(vExt)),vExt);
        hold on
    end   
    clearvars -except rate files f coefvar acStr SpkISI
 
end



% misc functions

% save variables as mouse number
% [~,name]=fileparts(files(1).folder); 
% save(name,'coefVar','rate','acStr') 

% used to combine matlab files into one file/treatment type
% clear
% rateAll=[];
% files=dir(fullfile('C:\Users\...','*.mat'));
% for i=1:length(files)
%     load(files(i).name)
%     rateAll=[rateAll rate];
%     clearvars -except files rateAll
% end
% [~,name]=fileparts(files(1).folder);
% save rateAll 'rateAll' 

% clear
% files=dir(fullfile('C:\Users\...','*.mat'));
% for i=1:length(files)
%     load(files(i).name)
%     errorbar(i+11,mean(rate),std(rate)./sqrt(size(rate,2)),'o','MarkerSize',10,'color','b');
%     clearvars -except files
% % end
 
% clear
% files=dir(fullfile('C:\Users\...','*.mat'));
% for i=1:length(files)
%     load(files(i).name)
%     rateAll(i)=mean(rate)
%     clearvars -except files rateAll
%     hold on
% end
