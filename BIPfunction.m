%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicubic interpolation and extrapolation iteration method for high %
%resolution digital holographic reconstruction

%Zhengzhong Huang, Liangcai Cao*
%State Key Laboratory of Precision Measurement Technology and Instruments, %Department of Precision Instruments, Tsinghua University, Beijing 100084, %China
%hzz19@mails.tsinghua.edu.cn
%clc@tsinghua.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang 2019
% The software is Matlab


function [out] = BIPfunction(in,num)

image=in;
[m,n]=size(in);
image=padarray(image,[2,2]);
image(1,:)=image(3,:);
image(2,:)=image(3,:);
image(n+3,:)=image(n+2,:);
image(n+4,:)=image(n+2,:);
image(:,1)=image(:,3);
image(:,2)=image(:,3);
image(:,n+3)=image(:,n+2);
image(:,n+4)=image(:,n+2);%padding image to [m+4,n+4]

for i=1:num*m
	u=rem(i,num)/num;%Get the corresponding decimal coordinates
	i1=floor(i/num)+2;
	A=[Weightfunction(1+u) Weightfunction(u) Weightfunction(1-u) Weightfunction(2-u)];   %the weight factor of four horizontal coordinates
	for j=1:num*n
        v=rem(j,num)/num;j1=floor(j/num)+2;%Get the corresponding decimal coordinates
        C=[Weightfunction(1+v);Weightfunction(v);Weightfunction(1-v);Weightfunction(2-v)];
        B=[image(i1-1,j1-1) image(i1-1,j1) image(i1-1,j1+1) image(i1-1,j1+2);...
           image(i1,j1-1) image(i1,j1) image(i1,j1+1) image(i1,j1+2);...
           image(i1+1,j1-1) image(i1+1,j1) image(i1+1,j1+1) image(i1+1,j1+2);...
           image(i1+2,j1-1) image(i1+2,j1) image(i1+2,j1+1) image(i1+2,j1+2);];
        out(i,j)=(A*B*C);
    end
end


