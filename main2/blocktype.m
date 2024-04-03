function out=blocktype(Tx_eff,Rx_eff,tab)
path = sprintf('%spilotParameters_Block_%s.mat','dataTx\',tab);
load(path);
H_pilot= Rx_eff(:,pilotSymbolIndx)./Tx_eff(:,pilotSymbolIndx);% Here is corresponding to block-type pilot,
                                                        %H_pilot is the channel transfer function of the pilot carrier.
PilotNumber=numel(pilotSymbolIndx);
H=zeros(size(Rx_eff,1),pilotSymbolIndx(end));
H(:,pilotSymbolIndx)=H_pilot;
for m=1:PilotNumber-1
    PilotInterval=pilotSymbolIndx(m+1)-pilotSymbolIndx(m);
    HcalculationIndex=pilotSymbolIndx(m);
    for l=1:PilotInterval-1
          H(:,HcalculationIndex+l)=(1-l/PilotInterval)*H_pilot(:,m)+(l/PilotInterval)*H_pilot(:,m+1);% Linear interpolationP
%         H(:,(L-1)*PilotInterval+1+l)=(H_pilot(:,L)+H_pilot(:,L+1))/2;% NO Linear interpolation
    end
end
CR=mean(H,2);% time domain mean is to reduce the amplitude noise?
Rx_eff=Rx_eff.*repmat(conj(CR)./(abs(CR).^2),1,size(Rx_eff,2));
% Rx_eff=Rx_eff.*(conj(H)./(abs(H).^2));
out=Rx_eff;
figure('name','Normalized response')
plot(20*log10(abs(CR)/abs(max(CR))),'-ro');ylabel('Normalized Roll-off[dB]');xlabel('SC index');
