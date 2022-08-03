clear; clc; close all;
addpath(genpath('C:\Users\asus\Desktop\tubitak_project'));
[SensEAT,path] = uigetfile('*.*','.mat');
SensEATFiles=load('-mat',SensEAT);
swallowData=SensEATFiles.datas;
breathingData=SensEATFiles.datas2;
ECGData=SensEATFiles.datas3;
ConsumeData=SensEATFiles.datas4;
%%  Swallow 
swallowData_DCoffset=swallowData-mean(swallowData);
swallowData_normalized=swallowData_DCoffset/max(swallowData_DCoffset);

figure(1)
plot(swallowData_normalized)
xlabel('Time')
ylabel('Amplitude')
title('Normalized Swallow Data')
% xlim([0 length(datas)])
% ylim([-1 1.1])
% [c,l]=wavedec(datas3,6,'sym4');
% a8 = wrcoef('a',c,l,'sym4',6);
% block_duration_datacorrected=datas3-a8;
fs=400;
% faxis=linspace(-fs/2,fs/2,length(datas3));
% plot(faxis,fftshift(abs(fft(datas3))));
block_duration_datacorrected=swallowData_normalized;
% Band    = (2 / fs) * [25, 199];
% [B, A]  = butter(6, Band, 'Bandpass');   
% fSignal = filtfilt(B, A, double(block_duration_datacorrected));
    % figure
Band    = (2 / fs) * [45, 55];
[B, A]  = butter(4, Band, 'stop');   
block_duration_datacorrected_last = filtfilt(B, A, double(block_duration_datacorrected));
   
N=length(block_duration_datacorrected_last);         %number of points
t=(0:N-1)/fs;   %time vector
sgf = sgolayfilt(block_duration_datacorrected_last,3,201);
% figure
% plot(sgf,'r')
% title(['Savitzky-Golay Smoothed - block: ' num2str(blocks{block_num})])
%     sgf=sgf(101:end-101);
ynew=sgf/max(sgf); 
% initialize filtered signal
eogF = ynew;
figure(2)
plot(eogF)
xlabel('Time')
ylabel('Amplitude')
%     xlim([1020.231193413892 260567.0467979074])
%     ylim([-0.0001904646550288437 0.0002626080637102741])
    % TKEO basic  % Teager–Kaiser energy operator to obtain EMG Bursts
for i=2:length(eogF)-1
    eogF(i) = ynew(i)^2 - ynew(i-1)*ynew(i+1);
end   
    % eogF=eogF(1:end-1);
[c,l] = wavedec(eogF,4,'sym4'); % 8 level decomposition  

    % for t=1:8
    %     D(:,t)=wrcoef('d',c,l,'db6',t);
    % end

    % low frequency components to get swallow patterns, so that approximation
    % of wavelets are obtained
clear A;
for t=1:4
    A(:,t)=wrcoef('a',c,l,'sym4',t);
end

A8=A(:,4);  % A8 is the filtered and swallow pattern obtained signal


A8_=A8(50:end-50);
figure(3)
plot(A8_)
xlabel('Time')
ylabel('Amplitude')
        rmsSwallows = sqrt(movmean(A8_.^2, 400));   % Burst Detection using RMS Value Over ‘WinLen’ Samples
figure(4)
plot(rmsSwallows)
title('Choose Threshold')
button = 1;
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
thresholdValue=(yy);  
close(figure(4)) 

x=400
for k=x:(length(rmsSwallows)-x)
    interval_=rmsSwallows(k-x+1:k+x);
    if(interval_(ceil(length(interval_)/2))==max(interval_) & interval_(ceil(length(interval_)/2))>thresholdValue)
        swallowpeaks(k) = rmsSwallows(k);
    end
end
swallowpeaks_pos=find(swallowpeaks>thresholdValue);
hv=ones(1,length(rmsSwallows));
swallowoffset=sqrt(-1)*hv;
swallowonset=swallowoffset;
    
for k=1:length(swallowpeaks_pos)
    mt=rmsSwallows(swallowpeaks_pos(k):swallowpeaks_pos(k)+450);
    mn=min(mt);
    swallowoffset(find(mt==mn)+swallowpeaks_pos(k)-1)=mn;
end
swallowoffset_last=find(swallowoffset~=sqrt(-1));

for m=1:length(swallowpeaks_pos)
    mt2=rmsSwallows(swallowpeaks_pos(m)-500:swallowpeaks_pos(m));
    mn2=min(mt2);
    swallowonset(swallowpeaks_pos(m)-500-1+find(mt2==mn2))= mn2;
end
swallowonset_last=find(swallowonset~=sqrt(-1));
    
figure(5)
plot(rmsSwallows)
hold on
plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
hold on
plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
hold on
plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
xlabel('Time')
ylabel('Amplitude')
title('Spontaneous Swallow Start')
    
SwallowDurations=swallowoffset_last-swallowonset_last;
for i=1:length(swallowpeaks_pos)
    triangle_(i,:)=[rmsSwallows(swallowonset_last(i)) rmsSwallows(swallowpeaks_pos(i)) rmsSwallows(swallowoffset_last(i))];
    areaUnderCurveFinal(i)=trapz(triangle_(i,:)); %trapz is best option to integrate 
end
SwallowAmplitudes=rmsSwallows(swallowpeaks_pos);
    
button = 1;
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
SpontaneousSwallowStart=(xx);  
close(figure(5)) 
figure(5)
plot(rmsSwallows)
hold on
plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
hold on
plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
hold on
plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
xlabel('Time')
ylabel('Amplitude')
title('Spontaneuos Swallow End')
button = 1;
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
SpontaneousSwallowEnd=(xx);
close(figure(5)) 

format longG
SpontaneousSwallowInterval=SpontaneousSwallowStart:SpontaneousSwallowEnd;
SpontaneousSwallowInterval=floor(SpontaneousSwallowInterval);
[tf,idx] = ismember(swallowpeaks_pos,SpontaneousSwallowInterval)
SpontaneousSwallows=find(idx>0);
SpontaneousSwallowRate=(SpontaneousSwallowEnd-SpontaneousSwallowStart)/length(SpontaneousSwallows);
SpontaneousSwallowsAmplitudes=SwallowAmplitudes(SpontaneousSwallows);
SpontaneousSwallowsDurations=SwallowDurations(SpontaneousSwallows);
SpontaneousSwallowsAreas=areaUnderCurveFinal(SpontaneousSwallows);

% Cued swallows
figure(5)
plot(rmsSwallows)
hold on
plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
hold on
plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
hold on
plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
xlabel('Time')
ylabel('Amplitude')
button = 1;
title('Cued Swallows Start')
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
CuedSwallowsStart=(xx);
close(figure(5)) 
figure(5)
plot(rmsSwallows)
hold on
plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
hold on
plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
hold on
plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
xlabel('Time')
ylabel('Amplitude')
button = 1;
title('Cued Swallows End')
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
CuedSwallowsEnd=(xx);
close(figure(5)) 

format longG
CuedSwallowInterval=CuedSwallowsStart:CuedSwallowsEnd;
CuedSwallowInterval=floor(CuedSwallowInterval);
[tf1,idx1] = ismember(swallowpeaks_pos,CuedSwallowInterval)
CuedSwallows=find(idx1>0);

CuedSwallowsAmplitudes=SwallowAmplitudes(CuedSwallows);
CuedSwallowsDurations=SwallowDurations(CuedSwallows);
CuedSwallowsAreas=areaUnderCurveFinal(CuedSwallows);

basefilename=[SensEAT(1:end-4) '_Swallows2'] ;
dlmwrite([basefilename '.csv'] ,[mean(CuedSwallowsAmplitudes)' mean(CuedSwallowsDurations/400)' mean(CuedSwallowsAreas)' SpontaneousSwallowRate' mean(SpontaneousSwallowsAmplitudes)' mean(SpontaneousSwallowsDurations/400)' mean(SpontaneousSwallowsAreas)'],'precision','%20.5f');
% Cued: 1 amplitude, 2 durations, 3 Area, 4 Spontaneous Rate, 5, Amplitude,
% 6, Duration, 7 Area
%%  Respiratory
fs=400
figure
plot(breathingData)
hold on
respiratoryRecording=breathingData-mean(breathingData);
plot(respiratoryRecording,'r')
[c,l]=wavedec(respiratoryRecording,10,'db6');
a16 = wrcoef('a',c,l,'db6',10);
respiratoryRecording_=respiratoryRecording-a16;

figure
plot(respiratoryRecording_)
% hold on
% plot(respiratoryData.resp,'r')

% plot(faxes,fftshift(abs(fft(respiratoryRecording_))));

fl=0.1; fu=0.5; % lower and upper cutoff freqs
wcl=2*fl/fs;
wcu=2*fu/fs;
N=100;
wn=[wcl wcu];
b=fir1(N,wn,'bandpass',hann(N+1)); %notch filter with hamming window
[h,w]= freqz(b,1,256);
j=1:length(respiratoryRecording_);
iv=zeros(1,N); %initialization vector for all taps to zero
respiratoryRecording_1=filter(b,1,respiratoryRecording_,iv);

 
respiratoryRecording_2 = smoothdata(respiratoryRecording_1,'SmoothingFactor',0.01);
figure(5);
plot(respiratoryRecording_2)
button = 1;
title('Choose Threshold')
while sum(button) <=1   % read input with right click only for threshold
    [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
end  % you can change
ThresholdForResp=(yy);
close(figure(5)) 
peaks=respiratoryRecording_2;

exhale=0;
x=100;

for k=x:(length(peaks)-x)
    gecici=peaks(k-x+1:k+x);
    if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>ThresholdForResp)
        exhale(k) = peaks(k);
    end
end
exhale_pos=find(exhale>0);
% for j=1:length(exhale_pos)-1
%     if(abs(exhale_pos(j)-exhale_pos(j+1))<=1000) 
%         if(exhale(exhale_pos(j))>exhale(exhale_pos(j+1)))
%             exhale(exhale_pos(j+1))=0;
%         elseif(exhale(exhale_pos(j))<exhale(exhale_pos(j+1)))
%             exhale(exhale_pos(j))=0;
% %             beatcount1=beatcount1+1;
%         end
%     end
%         
% end
% t=find(exhale>0);
% exhale_pos_last=t;

figure
plot(respiratoryRecording_2)
hold on
% plot(exhale_pos_last,respiratoryRecording_2(exhale_pos_last),'r+','MarkerFaceColor','r')
plot(exhale_pos,respiratoryRecording_2(exhale_pos),'r+','MarkerFaceColor','r')

hv=ones(1,length(respiratoryRecording_2));
inhale=sqrt(-1)*hv;
inhale2=sqrt(-1)*hv;

sbeatcount=0;
for k=1:length(exhale_pos)
    mt=respiratoryRecording_2(exhale_pos(k):exhale_pos(k)+750);
    mn=min(mt);
    inhale(find(mt==mn)+exhale_pos(k)-1)=mn;
%     sbeatcount=sbeatcount+1;
end
inhale_pos=find(inhale~=sqrt(-1));
% inhale_last=hv*sqrt(-1);

% for k=2:length(exhale_pos_last)
%     mt=respiratoryRecording_2(exhale_pos_last(k)-500:exhale_pos_last(k));
%     mn=min(mt);
%     inhale2(find(mt==mn)+exhale_pos_last(k)-501)=mn;
% %     sbeatcount=sbeatcount+1;
% end
% inhale_pos2=find(inhale2~=sqrt(-1));
% inhale_last2=hv*sqrt(-1);
% inhale_pos3=sort([inhale_pos,inhale_pos2]);
% 
% for i=1:length(inhale_pos3)-1
%     if(inhale_pos3(i+1)-inhale_pos3(i)<500)
%         inhale_pos3(i+1)=0;
%     end
% end
% indexofspos3=find(inhale_pos3>0);
% inhale_pos3=inhale_pos3(indexofspos3(1:end));
% 
% for k=1:length(inhale_pos3)
%     while(respiratoryRecording_2(inhale_pos3(k))>respiratoryRecording_2(inhale_pos3(k)+1))
%         inhale_pos3(k)=inhale_pos3(k)+1;
%     end
% end


hold on
plot(inhale_pos,respiratoryRecording_2(inhale_pos),'r*','MarkerFaceColor','r')

legend('Filtered data','Exhale Peaks','Inhale Troughs');

ylabel({'Normalized Nasal air-flow (°C)'});
xlabel({'Time (seconds)'})
format longG

BreathingAmplitude=respiratoryRecording_2(exhale_pos)-respiratoryRecording_2(inhale_pos);
BreathingDuration=inhale_pos-exhale_pos; % between inhale and exhale.
BreathingInterval=inhale_pos(2:end)-inhale_pos(1:end-1);
BreathingFreq=60./(BreathingInterval/400);

basefilename=[SensEAT(1:end-4) '_Respiratory'] ;
dlmwrite([basefilename '.csv'] ,[mean(BreathingAmplitude)' mean(BreathingDuration/400)' mean(BreathingInterval/400)' mean(BreathingFreq)'],'precision','%20.5f');
%1-Breathing Amplitude, 2-Breathing Duration between exhale-inhale,
%3-Breathing interval, between inhale to inhale, 4-Breathing Freq.


% meanofExhales=mean(respiratoryRecording_2(exhale_pos));
% apneas=find(respiratoryRecording_2(exhale_pos)<meanofExhales);
% exhale_pos(apneas)=[];
%% ECG ANALYSIS %%%%%%%%%%%%%%  ECGData

ECG_data=ECGData-mean(ECGData);
[c,l]=wavedec(ECG_data,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
ecgcorrected=ECG_data-a8;
%% Adaptive Filtering
fl=20;  % lower and upper cutoff freqs
wcl=2*fl/fs;
% faxes=linspace(-fs/2,fs/2,length(ecgsignal));
% plot(faxes,fftshift(abs(fft(ecgsignal))));
% axis([0 100 0 4e6])
N=50;
wn=[wcl];
b=fir1(N,wn,hamming(N+1)); %notch filter with hamming window
% [h,w]= freqz(b,1,256);
% j=1:length(ecgcorrected);
iv=zeros(1,N); %initialization vector for all taps to zero
ecgfiltered_last=filter(b,1,ecgcorrected,iv);
ecgfiltered_last = sgolayfilt(ecgfiltered_last,3,11);
% figure; plot(ecgsignal/max(ecgsignal)); title('Raw signal');
% figure; plot(ecgfiltered_last/max(ecgfiltered_last)); title('Filtered signal');
ecgfiltered_last=ecgfiltered_last/max(ecgfiltered_last);
figure
plot(ecgfiltered_last)
%% Wavelet Transform
% Using Wavelet Transform to Decompose signal
[c,l] = wavedec(ecgfiltered_last,8,'db6');
%% Wavelet Reconstruction with Coefficients
% DESCRIPTIVE TEXT
for t=1:8
    D(:,t)=wrcoef('d',c,l,'db6',t);
end
a2=wrcoef('a',c,l,'db6',2);
%% R-Peak Detection
% Use selected coefficients to R Peak Detection
e1=D(:,3)+D(:,4)+D(:,5);
e2=(D(:,4).*(D(:,3)+D(:,5)))/2.^8;
R_Peak_Detect_Ecg=e1.*e2;
R_Peak_Detect_Ecg_Positive = zeros(1,length(R_Peak_Detect_Ecg));
for k=1:length(R_Peak_Detect_Ecg)   
if R_Peak_Detect_Ecg(k)>0
R_Peak_Detect_Ecg_Positive(k)=R_Peak_Detect_Ecg(k);
end
end
last_ecg=R_Peak_Detect_Ecg_Positive;
threshold=max(last_ecg);
threshold=threshold*0.01;
r_peak=0;
x=ceil(fs/2);

for k=x:(length(last_ecg)-x)
        gecici=last_ecg(k-x+1:k+x);
        if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>threshold)
        r_peak(k) = last_ecg(k);
    end
end
r_peak_pos=find(r_peak>0);
% for j=1:length(r_peak_pos)-1
%     if(abs(r_peak_pos(j)-r_peak_pos(j+1))<=20) 
%         if(r_peak(r_peak_pos(j))>r_peak(r_peak_pos(j+1)))
%         r_peak(r_peak_pos(j+1))=0;
%         elseif(r_peak(r_peak_pos(j))<r_peak(r_peak_pos(j+1)))
%             r_peak(r_peak_pos(j))=0;
%         end
%     end
%         
% end
% t=find(r_peak>0);
% r_peak_pos=t;

r_peak_last=zeros(1,length(last_ecg));
for t=1:length(r_peak_pos)
    mt3=ecgfiltered_last(r_peak_pos(t)-23:r_peak_pos(t)+23);
    mn3=max(mt3);
    r_peak_last(find(mt3==mn3)+r_peak_pos(t)-24)=mn3;
end
r_peak_pos_last=find(r_peak_last>0);

% for i=1:length(r_peak_pos_last)-1
% if((r_peak_pos_last(i+1)-r_peak_pos_last(i))<150) & (ecgfiltered_last(r_peak_pos_last(i+1))>ecgfiltered_last(r_peak_pos_last(i)))
%         r_peak_pos_last(i)= [];
%     else if(((r_peak_pos_last(i+1)-r_peak_pos_last(i))<150) & (ecgfiltered_last(r_peak_pos_last(i+1))<ecgfiltered_last(r_peak_pos_last(i))))
%                 r_peak_pos_last(i+1)= [];
%         end
%     end
% end
%% Q & S Detection
hv=ones(1,length(ecgfiltered_last));
s_peak=sqrt(-1)*hv;
q_peak=s_peak;

for k=1:length(r_peak_pos_last)
    mt=ecgfiltered_last(r_peak_pos_last(k):r_peak_pos_last(k)+25);
    mn=min(mt);
    s_peak(find(mt==mn)+r_peak_pos_last(k)-1)=mn;
end
s_peak_pos_last=find(s_peak~=sqrt(-1));

for m=1:length(r_peak_pos_last)
    mt2=ecgfiltered_last(r_peak_pos_last(m)-25:r_peak_pos_last(m));
    mn2=min(mt2);
    q_peak(r_peak_pos_last(m)-25-1+find(mt2==mn2))= mn2;
end
q_peak_pos_last=find(q_peak~=sqrt(-1));

%% P and T Detection   ******* Burda kald?m Buraya kadar Süper 
e4=D(:,4)+D(:,5)+D(:,6)+D(:,7)+D(:,8);

t_peak=hv*sqrt(-1);
for t=1:length(s_peak_pos_last)
    mt6=e4(s_peak_pos_last(t):s_peak_pos_last(t)+80);
    mn6=max(mt6);
    t_peak(find(mt6==mn6)+s_peak_pos_last(t)-1)=mn6;
end
t_peak_pos=find(t_peak~=sqrt(-1));

t_peak_last=hv*sqrt(-1);

for t=1:length(t_peak_pos)
    mt9=ecgfiltered_last(t_peak_pos(t)-3:t_peak_pos(t)+3);
    mn9=max(mt9);
    in = find(mt9==mn9);
    t_peak_last(in(1)+t_peak_pos(t)-4)=mn9;
end

t_peak_pos_last = find(t_peak_last~=sqrt(-1));

% for x=1:length(t_peak_pos_last)
% if(r_peak_pos_last(x+1)-t_peak_pos_last(x)<(ceil(fs/bpm)*6))
%     t_peak_last(t_peak_pos_last(x))=0;
% mt13=ecgfiltered_last(s_peak_pos_last(x):s_peak_pos_last(x)+(ceil(fs/bpm)*8));
% mn13=max(mt13);
% in=find(mt13==mn13);
% t_peak_last(in(1)+s_peak_pos_last(x)-1)=mn13;
% end
% end
% t_peak_pos_last = find(t_peak_last~=sqrt(-1));

p_peak=hv*sqrt(-1);

for p=2:length(q_peak_pos_last)
    mt10=e4(q_peak_pos_last(p)- 45:q_peak_pos_last(p));
    mn10=max(mt10);
    p_peak(find(mt10==mn10)+q_peak_pos_last(p)-45-1)=mn10;
end

p_peak_pos=find(p_peak~=sqrt(-1));

p_peak_last=hv*sqrt(-1);

for p=1:length(p_peak_pos)
    mt11=ecgfiltered_last(p_peak_pos(p)- 3:p_peak_pos(p)+3);
    mn11=max(mt11);
    in = find(mt11==mn11);
    p_peak_last(in(1)+p_peak_pos(p)-3+1)=mn11;
end

p_peak_pos_last=find(p_peak_last~=sqrt(-1));
%%
%%P ve T Start and Final Points
% P Start
p_peak_start=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt15=ecgfiltered_last(p_peak_pos_last(i)-30:p_peak_pos_last(i));
    mn15=min(mt15);
    in = find(mt15==mn15);
    p_peak_start(in(1)+p_peak_pos_last(i)-30-1)=mn15;
end
p_peak_start_pos=find(p_peak_start~=sqrt(-1));

%P Final
p_peak_final=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt16=ecgfiltered_last(p_peak_pos_last(i):p_peak_pos_last(i)+15);
    mn16=min(mt16);
    in = find(mt16==mn16);
    p_peak_final(in(1)+p_peak_pos_last(i)-1)=mn16;
end
p_peak_final_pos=find(p_peak_final~=sqrt(-1));

% for i=1:length(p_peak_pos_last)-1
%     if(q_peak_pos_last(i+1)<=p_peak_final_pos(i))
%         q_peak_pos_last(i+1)=q_peak_pos_last(i+1)+3;
%         p_peak_final_pos(i)=p_peak_final_pos(i)-3;
%     end
% end

% T Final
t_peak_final=hv*sqrt(-1);
for i=1:length(t_peak_pos_last)
    mt18=ecgfiltered_last(t_peak_pos_last(i):t_peak_pos_last(i)+37);
    mn18=min(mt18);
    in = find(mt18==mn18);
    t_peak_final(in(1)+t_peak_pos_last(i)-1)=mn18;
end
t_peak_final_pos=find(t_peak_final~=sqrt(-1));

    t_peak_start=hv*sqrt(-1);   
for i=1:length(t_peak_pos_last)
    mt17=ecgfiltered_last(t_peak_pos_last(i)-45:t_peak_pos_last(i));
    mn17=min(mt17);
    in = find(mt17==mn17);
    t_peak_start(in(1)+t_peak_pos_last(i)-46)=mn17;
end
t_peak_start_pos=find(t_peak_start~=sqrt(-1));     

figure;
hold on;
plot(ecgfiltered_last)
plot(r_peak_pos_last,ecgfiltered_last(r_peak_pos_last),'r+','MarkerFaceColor','r')
plot(s_peak_pos_last,ecgfiltered_last(s_peak_pos_last),'r*','MarkerFaceColor','r')
plot(q_peak_pos_last,ecgfiltered_last(q_peak_pos_last),'r.','MarkerFaceColor','r')
plot(p_peak_start_pos,ecgfiltered_last(p_peak_start_pos),'k+','MarkerFaceColor','r')
plot(p_peak_pos_last,ecgfiltered_last(p_peak_pos_last),'k*','MarkerFaceColor','r')
plot(p_peak_final_pos,ecgfiltered_last(p_peak_final_pos),'k.','MarkerFaceColor','r')
plot(t_peak_start_pos,ecgfiltered_last(t_peak_start_pos),'b+','MarkerFaceColor','r')
plot(t_peak_pos_last,ecgfiltered_last(t_peak_pos_last),'b*','MarkerFaceColor','r')
plot(t_peak_final_pos,ecgfiltered_last(t_peak_final_pos),'b.','MarkerFaceColor','r')
% axis([0 5000 -0.4 1])
% set(gca,'FontName','Times New Roman','XTick',...
%     [0 1000 2000 3000 4000 5000 ],'XTickLabel',...
%     {'0','2.5','5','7.5','10','12.5'});
ylabel({'ECG signal (mV)'});
xlabel({'Time (seconds)'})
% figure
% plot(ECGData/max(ECGData))
% axis([0 5000 0.55 1])
% set(gca,'FontName','Times New Roman','XTick',...
%     [0 1000 2000 3000 4000 5000 ],'XTickLabel',...
%     {'0','2.5','5','7.5','10','12.5'});
% ylabel({'ECG signal (mV)'});
% xlabel({'Time (seconds)'})
HRVres=0;
for i=1:length(r_peak_pos_last)-1
    HRVres(i)=r_peak_pos_last(i+1)-r_peak_pos_last(i);
end
% HRV is calculated as the standard deviation of all of the 
% RR intervals (the distance between each “R” of the QRS complex).
hrv_last=std(HRVres);
HR=60*400*length(r_peak_pos_last)/length(ECGData);

st_interval= zeros(1,length(t_peak_start_pos));
for i =1:length(st_interval)-1
    st_interval(i)=t_peak_final_pos(i)-s_peak_pos_last(i); 
    %t bitiş - s
end
mean_st_interval= mean(st_interval)/400;
st_segment=zeros(1,length(t_peak_pos_last));
for i =1:length(st_interval)
    st_segment(i)=t_peak_start_pos(i)-s_peak_pos_last(i);
    %t başlangıç - s
end
mean_st_segment=mean(st_segment)/400;
% QRS Interval and amplitude
qrs_amp=zeros(1,length(r_peak_pos_last));
qrs_int=zeros(1,length(q_peak_pos_last));
for i=1:length(qrs_int)
    qrs_int(i)=s_peak_pos_last(i)-q_peak_pos_last(i);
end
for i=1:length(qrs_amp)
    qrs_amp(i)=ecgfiltered_last(r_peak_pos_last(i))-ecgfiltered_last(s_peak_pos_last(i));
end
mean_qrs_int=mean(qrs_int)/400;
mean_qrs_amp=mean(qrs_amp);
rr_int=zeros(1,length(r_peak_pos_last)-1);

for i=1:length(r_peak_pos_last)-1
    rr_int(i)=(r_peak_pos_last(i+1)-r_peak_pos_last(i));
end
rr_intmean=mean(rr_int)/400;
% QT interval 
qt_int=zeros(1,length(q_peak_pos_last));

for i=1:length(q_peak_pos_last)-2
    qt_int(i)=t_peak_final_pos(i)-q_peak_pos_last(i);
end

qt_intmean=mean(qt_int)/400;

Qtc_last=0;
for i=1:length(r_peak_pos_last)-1
    Qtc_last(i)=(qt_int(i)/rr_int(i))^(1/3);
end
Qtc_mean=mean(Qtc_last);

basefilename=[SensEAT(1:end-4) '_ECG'] ;
dlmwrite([basefilename '.csv'] ,[hrv_last' HR' mean_st_interval' mean_st_segment' mean_qrs_int' mean_qrs_amp' rr_intmean' qt_intmean' Qtc_mean'],'precision','%20.5f');

%1-HRV, 2-HR, 3-ST interval, 4- ST segment, 5- QRS int, 6- QRS amp, 
%7- RR int, 8- QT int, 9- QTc
%% WEIGHT
%%%%%%%%%%%%% WEIGHT MEASUREMENT: ConsumeData %%%%%%%%%%%%%%%%%%% 
figure
plot(ConsumeData)
xlabel('Time (Seconds)')
ylabel('Weight (grams)')
CuedConsumation=ConsumeData(floor(CuedSwallowsStart))-ConsumeData(floor(CuedSwallowsEnd));
SpontaneousConsumation=ConsumeData(floor(SpontaneousSwallowStart))-ConsumeData(floor(SpontaneousSwallowEnd));
ti_CC=floor(CuedSwallowsEnd-CuedSwallowsStart)/fs;
ti_SC=floor(SpontaneousSwallowEnd-SpontaneousSwallowStart)/fs;

SC_Rate=SpontaneousConsumation/ti_SC;

basefilename=[SensEAT(1:end-4) '_Consumption'] ;
dlmwrite([basefilename '.csv'] ,[CuedConsumation' SpontaneousConsumation' ti_CC' ti_SC' SC_Rate'],'precision','%20.5f');
% 1- CuedConsumationAmount 2- SpontaneousConsumationAmount
% 3-TimeIntervalCuedConsump, 4- TimeIntervalSpontCons 5- SpontConsRate
