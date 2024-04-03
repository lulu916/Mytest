H1=uint8(imread('lena2.bmp'));  
subplot(2,2,1);  
imshow(H1);  
title('原图');  
subplot(2,2,2);  
imhist(H1);  
title('原图直方图'); 

H2=adapthisteq(uint8(imread('jiamitu.bmp')));   
subplot(2,2,3);  
imshow(H2);  
title('加密图');  
subplot(2,2,4);  
imhist(H2);  
title('加密直方图'); 