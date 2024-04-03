function []=PAPRCCDF(TimeSigMatrix)
time_signal_matrix_in=(abs(TimeSigMatrix)).^2;
Max_time_signal_in=max(time_signal_matrix_in,[],1);
Mean_time_signal_in=mean(time_signal_matrix_in,1);
PAPR_Orignal=10*log10(Max_time_signal_in./Mean_time_signal_in);
zdBs=[4:0.1:12];
N_zdBs=length(zdBs);
for i=1:N_zdBs 
     CCDF_Orn(i)=sum(PAPR_Orignal>zdBs(i))/size(TimeSigMatrix,2);
end
figure('name','The PAPR CCDF Curve')
semilogy(zdBs(1:end),CCDF_Orn(1:end),'m:>')