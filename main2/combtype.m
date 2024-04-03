function out=combtype(Tx_eff,Rx_eff,Label)
path = sprintf('%sPilotParameters_Comb_%s.mat','dataTx\',Label);
load(path);

PilotHmatrix= Rx_eff(PilotCarrierIndex,:)./Tx_eff(PilotCarrierIndex,:);% Here is corresponding to block-type pilot,
                                                        %H_pilot is the channel transfer function of the pilot carrier.
PilotNumber=numel(PilotCarrierIndex);
H=zeros(PilotCarrierIndex(end),size(Rx_eff,2));
H(PilotCarrierIndex,:)=PilotHmatrix;
for m=1:PilotNumber-1
    pilotInterval=PilotCarrierIndex(m+1)-PilotCarrierIndex(m);
    for l=1:pilotInterval-1
        H(PilotCarrierIndex(m)+l,:)=(1-l/pilotInterval)*PilotHmatrix(m,:)+(l/pilotInterval)*PilotHmatrix(m+1,:);
    end
end
CR=mean(H,2);%Time domain mean is to reduce amplitude noise?
Rx_eff=Rx_eff.*repmat(conj(CR)./(abs(CR).^2),1,size(Rx_eff,2));
out=Rx_eff;
figure('name','Constellation')
plot(Rx_eff(:),'.')