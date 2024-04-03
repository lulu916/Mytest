%----------------------- no ofdm-------------------------
clc
clear all
P=imread('lena2.bmp');
P=double(P);
N=size(P,1);
CR=0.5;%—πÀı¬ 
m=CR*N;
lamuda=2.3;
%----------------- Chaos initial value----------------%
[x,y]=LSCM(0.99, 0.8, 0.5, 256);
u=abs((x+y)/2);
%u=u(1:256);
R=zeros(m,N);
R(1,1:N)=u;
for j=2:m
    R(j,:)=circshift(R(j-1,:),[0,1]);%—≠ª∑æÿ’Û
    R(j,1)=lamuda*R(j-1,N);
end
%---------------------CS process--------------------%
%ww=DWT(256);
%X1=ww*sparse(P)*ww';
%I1=full(X1);
figure(1),imshow(uint8(P));
I1=dct(P);  %DCT œ° ËªØ
%num = N*N;
%P1 = zigzagCoder(I1, num,lx0,ly0);
%I2=(reshape(P1,[],N))';
Y=R*I1; %—πÀı
YS=Y;
format long g
min=min(min(YS))
max=max(max(YS))
YSS=sum(sum(YS<0));
figure(2),imshow(uint8(Y));
imwrite(Y,'jiamitu.bmp');
Y(Y<0)=0;
%Y((Y>4096))=4096; 
Y=round(Y);

%---------------------ICS process----------------------%
R(1,1:N)=u;
for j=2:m
    R(j,:)=circshift(R(j-1,:),[0,1]);%—≠ª∑æÿ’Û
    R(j,1)=lamuda*R(j-1,N);
end
S=SL0(R,Y,N);
S1=idct(S);
%S1=ww'*sparse(S)*ww;
%S1=full(S1);
figure(3),imshow(uint8(S1));
imwrite(uint8(S1),'lena3.bmp');
%---------------------yuanshituxiang before CS----------------------%
IbitSequence1=reshape(P',[],1);
IbitSequence2=de2bi(IbitSequence1,'left-msb');
IbitSequence=reshape(IbitSequence2',1,[]);
