clc
clear all
P=imread('lena2.bmp');
P=double(P);
N=size(P,1);
CR=0.5;%压缩率
m=CR*N;
lamuda=2.3;
%----------------- Chaos initial value----------------%
[x,y]=LSCM(0.99, 0.8, 0.5, 256);
u=abs((x+y)/2);
%u=u(1:256);
R=zeros(m,N);
R(1,1:N)=u;
for j=2:m
    R(j,:)=circshift(R(j-1,:),[0,1]);%循环矩阵
    R(j,1)=lamuda*R(j-1,N);
end
%---------------------CS process--------------------%
%ww=DWT(256);
%X1=ww*sparse(P)*ww';
%I1=full(X1);
figure(1),imshow(uint8(P));
I1=dct(P);  %DCT 稀疏化
%num = N*N;
%P1 = zigzagCoder(I1, num,lx0,ly0);
%I2=(reshape(P1,[],N))';
Y=R*I1; %压缩
YS=Y;
YSS=sum(sum(YS<0));
figure(2),imshow(uint8(Y));
imwrite(Y,'jiamitu.bmp');
Y(Y<0)=0;
%Y((Y>4096))=4096; 
Y=round(Y);
bitSequence1=reshape(Y',[],1);
bitSequence2=de2bi(bitSequence1,'left-msb');
bitSequence=reshape(bitSequence2',1,[]);
save('tuxiangxulie.mat', 'bitSequence');
%---------------------OFDM  TX-------------------------%

%---------------------OFDM  RX-------------------------%

%---------------------ICS process----------------------%
load('Rtuxiangxulie.mat');
B1 = vec2mat(B3,13); 
B2 = (bi2de(B1,'left-msb'))'; %二进制转十进制
Y2=(reshape(B2,[],CR*N))';

R(1,1:N)=u;
for j=2:m
    R(j,:)=circshift(R(j-1,:),[0,1]);%循环矩阵
    R(j,1)=lamuda*R(j-1,N);
end
S=SL0(R,Y2,N);
S1=idct(S);
figure(3),imshow(uint8(S1));
imwrite(uint8(S1),'lena3.bmp');
%---------------------yuanshituxiang before CS----------------------%
IbitSequence1=reshape(P',[],1);
IbitSequence2=de2bi(IbitSequence1,'left-msb');
IbitSequence=reshape(IbitSequence2',1,[]);
img1 = double(imread('lena2.bmp'));
img2 = double(imread('lena3.bmp'));
%----------------------------mssim------------------------------------%
K = [0.01 0.03];
winsize = 11;
sigma = 1.5;
window = fspecial('gaussian', winsize, sigma);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'product';

msssim = ssim_mscale_new(img1, img2, K, window, level, weight, 'product');




