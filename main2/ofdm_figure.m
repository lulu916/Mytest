%==============================  Prameters  ===============================
clc
clear all
modOrder=16;% the modulation order
bitNmb=log2(modOrder);
datCarrierNmb=120;
fftSize=256;
CP=1/16;              %循环前缀长度占OFDM信号的1/16
awgSampleRate = 12.5*10^9;               % Gsa/s
dpoSampleRate =25*10^9;             % Gsa/s
P_mean_sqrt = 3.162277660168380;
pilotInsertType='no';
subShift=1;%subcarrier shift
carrierCut1=96;
carrierCut2=256;%from experiment
N=256;
% datCarrierIndx = [subShift+1:subShift+carrierCut1-1,subShift+carrierCut1+1:subShift+carrierCut2-1,subShift+carrierCut2+1:subShift+datCarrierNmb+2];
% datCarrierIndx = [subShift+1:subShift+carrierCut1-1,subShift+carrierCut1+1:subShift+datCarrierNmb+1];
datCarrierIndx = [subShift+1:subShift+datCarrierNmb];
%--------------------extract initial tuxiangxulie------------------------
load('tuxiangxulie.mat');
% try to reduce Peak-to-Average Ratio 24 512 1024,列移0行移24
%---------------   number of OFDM symbol   ---------------
symbolNmb =1000;% symbolNmb setup 较之前提升6倍以搭载Cs
if length(bitSequence) < datCarrierNmb*symbolNmb*bitNmb
   bitSequence = repmat(bitSequence,1,ceil(datCarrierNmb*symbolNmb*bitNmb/length(bitSequence)));
end
bitSequence = bitSequence(1,1:datCarrierNmb*symbolNmb*bitNmb);

%=========================  Modulation and S/P  ===========================
bitMatrix = vec2mat(bitSequence,bitNmb);          %向量转为矩阵 bitNmb=4 列
decSequence = bi2de(bitMatrix,'left-msb');        %二进制转十进制
modDat = qammod(decSequence,modOrder,[],'gray');  % 16 QAM 格雷星座映射
modDat=modDat/P_mean_sqrt;                        %归一化处理
datMatrix = vec2mat(modDat,datCarrierNmb).';
%---------------------- The data for BER calculation ----------------------
    berCalDat=datMatrix;
    path = sprintf('%sberCalDat.mat','dataTx\');
    save(path, 'berCalDat', '-v7.3'); % MAT 7.3 Can save > 2GB
% % %************************************************************************%
% %=====================    brownian  diffusion      =======================
%     plainDat=datMatrix;
%     diffusionDat=Diffusion(plainDat,bitNmb);
%     datMatrix=diffusionDat;
%=====================    brownian Motion  Encryption    ==================
    plainDat=datMatrix;
%     encDat=brownianMotionEnc(plainDat);
%     datMatrix=encDat;
% %-------------------------- Power normalized ----------------------------
%     datMatrix=datMatrix/P_mean_sqrt; 
%************************************************************************%        
if  strcmpi(pilotInsertType,'Comb')   %pilotInsertType=no
    pilotCarrierNmb=5;
    pilotCarrierIndx=ceil(linspace(1,numel(datCarrierNmb),pilotCarrierNmb));
    %pilot carrier index,This is because that the hermitian symmetric is performed.
    % 1 stands the DC carrier,which is zero padding.
    valDatCarrierIndx =1:numel(datCarrierNmb);
    valDatCarrierIndx(pilotCarrierIndx)=[];
    %-------------------------------CombType insertion---------------------
    % CombPilotSymbol=repmat(transmitSymbol(1,:),numel(PilotCarrierIndex),1);
    combPilotSymbol=datMatrix(pilotCarrierIndx,:);
    datPilotMatrix(pilotCarrierIndx,:)=combPilotSymbol;
    datPilotMatrix(valDatCarrierIndx,:)=datMatrix(valDatCarrierIndx,:);
    datPilotMatrix=datMatrix;
    path = sprintf('%spilotParameters_Comb.mat','dataTx\');
    save(path,'pilotCarrierIndx','valDatCarrierIndx');
elseif strcmpi(pilotInsertType,'Block')
    pilotSymbolNumber=19;
    pilotSymbolIndx=ceil(linspace(1 ,symbolNmb,pilotSymbolNumber));
    datSymbolIndx=1:symbolNmb;
    datSymbolIndx(pilotSymbolIndx)=[];
    %-------------------------------BlockType insertion------------------------
    blockPilotSymbol=blkPilotSymbGen2(datCarrierNmb,pilotSymbolNumber,13);
    datPilotMatrix(:,pilotSymbolIndx)=blockPilotSymbol;
    datPilotMatrix(:,datSymbolIndx)=datMatrix(:,datSymbolIndx);
    datPilotMatrix=datMatrix;
    path = sprintf('%spilotParameters_Block.mat','dataTx\');
    save(path,'pilotSymbolIndx','datSymbolIndx');
else
    datPilotMatrix=datMatrix;
end
%---------------------- The data for channel estimate ---------------------
    chnEstimateDat=datPilotMatrix;
    path = sprintf('%schnEstimateDat.mat','dataTx\');
    save(path, 'chnEstimateDat', '-v7.3'); % MAT 7.3 Can save > 2GB  
    
%===========================  OFDM modulation  ============================
transmissionDat=chnEstimateDat;
clear fftInMatrix
%---------------------------  Hermitian Symmetry  -------------------------
hermitIndex=zeros(fftSize/2-1,symbolNmb);   %fftSize=256;symbolNmb =200; 
hermitIndex(datCarrierIndx-1,:)=transmissionDat;

%-----------------  One Hermitian structure  -----------------
fftInMatrix=[zeros(1,symbolNmb);hermitIndex;zeros(1,symbolNmb);flipud(conj(hermitIndex))];%??/
%---------------------------         IFFT       ---------------------------
timDatMatrix= ifft(fftInMatrix); 
%---------------------------  PAPR CCDF curve  ----------------------------
PAPRCCDF(timDatMatrix)

%======================  Add Cyclic Prefix and P/S   ======================
timDatMatrixCP = [timDatMatrix(size(timDatMatrix,1)*(1-CP)+1:end,:);timDatMatrix];       % Add CP
timSignalSeries = reshape(timDatMatrixCP,1,size(timDatMatrixCP,1)*size(timDatMatrixCP,2)); % Change to row vector P\S
datLen=length(timSignalSeries);
timSignalSeries = timSignalSeries - mean(timSignalSeries); % Remove the DC 
%--------------------------------headpoint------------------------------------------
head       = load('PRBS7.txt');                                 % length=128
head       = head.';
syncHeader = kron(head,ones(1,2*2));                           % length=128*2*Nb
syncHeader = [syncHeader syncHeader] - mean(syncHeader);        % length=128*2*Nb*2 and make sure that mean(synchead)=0;
transmissionSignal    = [syncHeader*max(abs(timSignalSeries))/2 timSignalSeries];
headerLen=length(syncHeader);
%----------------------- Stend to AWG by mat file -------------------------
    awg_Tx = transmissionSignal.';
    x = length(awg_Tx); 
    t = 1:1:x; 
    baseMarkers = uint8(square(2*pi*1/x*t,50));
    baseWfm_I = awg_Tx; % Generate Sine Wave
    Waveform_Name_1 = ['OFDM']; 
    Waveform_Data_1 = baseWfm_I'; %already a double array 
    Waveform_M1_1 = baseMarkers; %already uint8 array 
    Waveform_M2_1 = baseMarkers;
    %save the output data
    path = sprintf('%sOFDM.mat','dataTx\');
    save(path, '*_1', '-v7.3'); % MAT 7.3 Can save > 2GB
    
%*******************解密解密解密解密解密解密解密解密解密**********************    
%*******************解密解密解密解密解密解密解密解密解密**********************
% path = sprintf('%sofdmParameter.mat','dataTx\');
% load(path);
awgSampleRate = awgSampleRate*10^9;            % Gsa/s
dpoSampleRate =dpoSampleRate*10^9;             % Gsa/s
% % pdInputPower=opticalPower ;
%--------------------------  Variable initialize  -------------------------
% % berTotalCase=zeros(1,sampleTimes);
% % berU1Case=zeros(1,sampleTimes);
% % berU2Case=zeros(1,sampleTimes);    

        path = sprintf('%sOFDM.mat','dataTx\');
        load(path);
        transmissionDat=repmat(Waveform_Data_1,1,4);
        transmissionDat=circshift(transmissionDat,[0 2000]);
        %transmissionDat=awgn(transmissionDat,18,'measured');%新加的,欲使函数将在加入噪声之前测定信号强度则加measured
        reSampleDat=transmissionDat;
        
 %===================  Find the Synchronization header  ====================
for i=1:length(reSampleDat)/2  %128 is synchronization header length
    correlationCoeficient(i)=sum((reSampleDat(i:i+headerLen/2-1)).*conj((reSampleDat(i+headerLen/2:i+headerLen-1))));
end
[peak,syncHeader]=max(abs(correlationCoeficient));
% % --------------------  Synchronization header plot  --------------------
%figure('name','Synchronization Header')
%plot(abs(correlationCoeficient));
%------------------  the synchronization output data  ---------------------
syncHeader=syncHeader+headerLen;%chang this parameter to obtain effective results,
Rx_series = reSampleDat(syncHeader:datLen+syncHeader-1);
% %========================   OFDM demodulation   =========================
CPshift = 2;
%-------------------------  S/P and Remove CP  ----------------------------
Rx_parallel = vec2mat(Rx_series,fftSize*(1+CP)).';
Rx_parallel = Rx_parallel(fftSize*CP+1-CPshift:fftSize*(1+CP)-CPshift,:); % Discard the CP
%---------------------------------  FFT  ----------------------------------
Rx_parallel = fft(Rx_parallel,fftSize)/(fftSize/2);
%--------------------  extracted the transmission data --------------------
Rx_parallel = Rx_parallel(datCarrierIndx,:);


%========================    Channel estimate   ===========================
path = sprintf('%schnEstimateDat.mat','dataTx\');
load(path); % MAT 7.3 Can save > 2GB
Tx_eff = chnEstimateDat;
Rx_eff = Rx_parallel;
if strcmpi(pilotInsertType,'Comb')
%-----------------------------Comb type pilot----------------------
    Rx_eff=combtype(Tx_eff,Rx_eff,tag);
elseif strcmpi(pilotInsertType,'Block')
%-----------------------------Block type pilot---------------------
    Rx_eff=blocktype(Tx_eff,Rx_eff,tag);
else
%------------------------Original compensation code--------------------
    npilot=Rx_eff./Tx_eff;
    CR=mean(npilot,2);
    Rx_eff=Rx_eff.*repmat(conj(CR)./(abs(CR).^2),1,symbolNmb);
end
%*************************************************
singalPower=mean(abs(Tx_eff.^2),2);
noisePower=var(Rx_eff-Tx_eff,1,2);
SNR=10*log10(singalPower./noisePower);
figure('name','The SNR distribution.')
plot(SNR)
SNR_total=mean(SNR);
% =====================    brownian Motion  decryption    ==================
    chiperDat=Rx_eff;
%    Rx_parallel=brownianMotionDec(chiperDat);
%    Rx_eff=Rx_parallel;
% % =====================     brownian  deDiffusion     ==================
% %     chiperDat=Rx_eff;
% %     chiperDat=qammod(chiperDat,2^bitNmb,[],'gray');
% %    Rx_eff=deDiffusion(chiperDat,bitNmb);
%     Rx_parallel=qamdemod(deDiffDat,2^bitNmb,[],'gray');
%      Rx_parallel=qamdemod(deDiffDat,2^bitNmb,'gray');   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================   QFactor calculate   ==========================
for mm = 1:length(datCarrierNmb)
    Q_factor(mm) = QFactor(Rx_eff(mm,:), Tx_eff(mm,:));
end
min_Qfactor=min(Q_factor);
min_nml=find(Q_factor==min_Qfactor);
% figure('name','The Q_factor distribution.')
% plot(Q_factor)
%===========================  total bit rate   ============================
Bitrat_CH1 = datCarrierNmb*bitNmb*awgSampleRate/fftSize/(1+CP);
%========================   Calculate the BER   ===========================
%--------------------------  the total BER  -------------------------------
path = sprintf('%sberCalDat.mat','dataTx\');
load(path); 
Tx_eff=berCalDat*P_mean_sqrt;
Rx_eff=Rx_eff*P_mean_sqrt;
Tx_eff_decode = qamdemod(Tx_eff,2^bitNmb,[],'gray');
Rx_eff_decode = qamdemod(Rx_eff,2^bitNmb,[],'gray');
%--------------------extract initial tuxiangxulie-------------------------
Rx_eff_demod=reshape(Rx_eff_decode,[],1);%变为1列
Rx_eff_bi=de2bi(Rx_eff_demod,'left-msb');
Rx_eff_bitSequence=reshape(Rx_eff_bi.',1,[]);
CR=0.5;
B3=Rx_eff_bitSequence(1,1:CR*N*N*13);
save('Rtuxiangxulie.mat', 'B3');
%------------------ recover  Rx_eff_decode--------------------------------
if length(B3) < datCarrierNmb*symbolNmb*bitNmb
  B7= repmat(B3,1,ceil(datCarrierNmb*symbolNmb*bitNmb/length(B3)));
end
B7 = B7(1,1:datCarrierNmb*symbolNmb*bitNmb);
B7Matrix = vec2mat(B7,bitNmb); 
B7decSequence = bi2de(B7Matrix,'left-msb');  
Rx_eff_decode = vec2mat(B7decSequence,datCarrierNmb).';
%------------------  Plot the Rx_eff constellation  ----------------------
figure('name','Constellation')
refpts = qammod((0:(2^bitNmb-1))',2^bitNmb,[],'gray');
colorVal=linspecer(2^bitNmb,'sequential');
for indx=1:2^bitNmb
    plot(Rx_eff(Rx_eff_decode==(indx-1)),'.')
    colormap Lines
    hold on
end
plot(refpts,'r*')   
%--------------------- Another demodulation method  ---------------------
[~,berTotal,EI] = biterr(Tx_eff_decode,Rx_eff_decode);
berTotalCase(t)=berTotal;
berPerCarrier=sum(EI,2)/(bitNmb*datCarrierNmb);
figure('name','BER distribution')
plot(datCarrierIndx,berPerCarrier)
%----------------------------   SNR-BER   ---------------------------------
%SNR = [7:1:19]; %  CR=0.25
%BER1 = [-0.91362 -0.9940 -1.0975 -1.2075 -1.3552 -1.5233 -1.7242 -1.9814 -2.2768 -2.6684 -3.2226 -3.7648 -4.2833 ];
%BER2=  [-0.92167 -1.0016 -1.0991 -1.2146 -1.3513 -1.5177 -1.7229 -2.0011 -2.2864 -2.7380 -3.1043 -3.7270 -4.1694 ];
%BER3= -0.001*randi([470 530],1,13);
SNR = [7:1:19]; %  CR=0.5
BER1 = [-0.92462 -1.0014 -1.1053 -1.2337 -1.3635 -1.5366 -1.7447 -1.9831 -2.3139 -2.6761 -3.2041 -3.8617 -4.3802 ];
BER2=  [-0.92876 -1.0077 -1.1012 -1.2215 -1.3626 -1.5276 -1.7378 -2.0035 -2.3072 -2.6856 -3.1298 -3.6726 -4.9031 ];
BER3= -0.001*randi([470 530],1,13);
%SNR = [7:1:19]; %  CR=0.75
%bitrates=2.2059e+19 ;
%BER1 = [-0.92506 -1.0135 -1.1096 -1.2358 -1.3716 -1.5462 -1.7365 -2.0070 -2.3288 -2.6951 -3.1389 -3.7134 -4.6812 ];
%BER2=  [-0.92729 -1.0103 -1.1119 -1.2253 -1.3652 -1.5476 -1.7611 -2.0127 -2.3206 -2.7114 -3.2125 -3.8542 -4.5721 ];
%BER3= -0.001*randi([470 530],1,13);
%SNR = [7:1:19];  %  CR=1
%BER1 = [-0.93459 -1.0180 -1.1140 -1.2337 -1.3786 -1.5548 -1.7652 -2.0239 -2.3506 -2.7627 -3.1707 -3.8819 -4.7270 ];
%BER2=  [-0.93413 -1.0165 -1.1184 -1.2215 -1.3832 -1.5551 -1.7696 -2.0359 -2.3512 -2.7452 -3.2850 -3.8300 -4.3802 ];
%BER3= -0.001*randi([470 530],1,13);
figure('name','ONU 在不同SNR下接收不同OFDM信号的BER曲线')
plot(SNR(1:end),BER1(1:end),'m:>');
hold on
plot(SNR(1:end),BER2(1:end),'k:o');
hold on
plot(SNR(1:end),BER3(1:end),'c:o');
axis([7 19 -5 0]);
xlabel('SNR') ;
ylabel('log10(BER)') ;
%----------------------------  User BER   ---------------------------------
dataCarrierU1=1:2:datCarrierNmb;
dataCarrierU2=2:2:datCarrierNmb;
[~,berU1,~] = biterr(Tx_eff_decode(dataCarrierU1,:),Rx_eff_decode(dataCarrierU1,:),bitNmb);
[~,berU2,~] = biterr(Tx_eff_decode(dataCarrierU2,:),Rx_eff_decode(dataCarrierU2,:),bitNmb);
berU1Case(t)=berU1;
berU2Case(t)=berU2;
format long g
% % ==================  The results display  =====================
disp('*******************************************************************');
disp(['The SNR is ',num2str(SNR_total,'%2.4f')]);
disp(['The total bitrates is ',num2str(Bitrat_CH1,'%.4e')]);
disp(['The total BER is ',num2str(mean(berTotalCase),'%.4e')]);
disp(['The ONU1 BER is ',num2str(log10(mean(berU1Case)),'%.4e')]);
disp(['The ONU2 BER is ',num2str(log10(mean(berU2Case)),'%.4e')]);
disp('*******************************************************************');
str1=['The SNR is ',num2str(SNR_total,'%2.4f')];
str2=['The total bitrates is ',num2str(Bitrat_CH1,'%.4e')];
str3=['The total BER is ',num2str(mean(berTotalCase),'%.4e')];
str4=['The ONU1 BER is ',num2str(mean(berU1Case),'%.4e')];
str5=['The ONU2 BER is ',num2str(mean(berU2Case),'%.4e')];
display=strvcat(str1,str2,str3,str4,str5);

