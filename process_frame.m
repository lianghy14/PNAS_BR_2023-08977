function [] = process_frame(dot,t,i_frame)
global v
%G0: 0,0; 600,0; 600,0; 600,600
%G1: 446,434; 820,434; 621,617; 1151,617
%G2: 

dir_fig = 'C:\Users\HLiang\Desktop\Case study LP\VideoFrames\';

f = read(v,i_frame);
% dot=[446,434; 621,617; 1151,617;820,434];       %取四个点，依次是左上，左下，右下,右上
y=[dot(1,1),dot(2,1),dot(3,1),dot(4,1)];        %四个原顶点
x=[dot(1,2),dot(2,2),dot(3,2),dot(4,2)];
h = 600; w = 600;
%这里是新的顶点，我取的矩形,也可以做成其他的形状
%大可以原图像是矩形，新图像是从dot中取得的点组成的任意四边形.:)
Y=[1,1,w,w];     
X=[1,h,h,1];

B=[X(1),Y(1),X(2),Y(2),X(3),Y(3),X(4),Y(4)]';   %变换后的四个顶点，方程右边的值
%联立解方程组，方程的系数
A=[x(1),y(1),1,0,0,0,-X(1)*x(1),-X(1)*y(1);             
   0,0,0,x(1),y(1),1,-Y(1)*x(1),-Y(1)*y(1);
   x(2),y(2),1,0,0,0,-X(2)*x(2),-X(2)*y(2);
   0,0,0,x(2),y(2),1,-Y(2)*x(2),-Y(2)*y(2);
   x(3),y(3),1,0,0,0,-X(3)*x(3),-X(3)*y(3);
   0,0 ,0,x(3),y(3),1,-Y(3)*x(3),-Y(3)*y(3);
   x(4),y(4),1,0,0,0,-X(4)*x(4),-X(4)*y(4);
   0,0,0,x(4),y(4),1,-Y(4)*x(4),-Y(4)*y(4)];%求解变换矩阵的行列式
cut_picture=zeros(h,w);
OutPicture=zeros(h,w,3);
OutPicture=cast(OutPicture,'uint8');%类型变换为uint8
fa=inv(A)*B;
fa=reshape([fa;1],[3,3]);%fa即是变换矩阵
cut=zeros(h+1,w+1);
BW=roipoly(f,y,x);%生成所选中区域的二值化图像cut
[coordinate_x,coordinate_y]=find(BW==1);%获得选中区域的坐标位置
coordinate=[coordinate_x,coordinate_y,ones(length(coordinate_y),1)];%增加z=1值以便进行左边变换
cut=coordinate*fa;%进行坐标变化
cut(:,1)=cut(:,1)./cut(:,3);
cut(:,2)=cut(:,2)./cut(:,3);%获得图像上对应的x，y值
cut=floor(cut+eps);
cut(find(cut(:,1:2)==0))=1;%边界为零的点变为1，
cut_coordinate=cut(:,1:2);%获得变换后得坐标位置
for k=1:3
    img=f(:,:,k);
    for i=1:length(coordinate)
        cut_picture(cut_coordinate(i,1),cut_coordinate(i,2))=img(coordinate(i,1),coordinate(i,2));
    end%讲对应坐标的像素对应
    cut_picture=cast(cut_picture,'uint8');
    B=zeros(5,5);
    D = padarray(cut_picture,[3 3],'replicate','both');%填充像素
    for i=4:h+3
        for j=4:w+3
            if D(i,j)==0
                no_zero=find(D(i-1:i+1,j-1:j+1)~=0);
                [m,n]=size(no_zero);
                if m~=0
                     B=D(i-1:i+1,j-1:j+1);
                else
                    no_zero=find(D(i-2:i+2,j-2:j+2)~=0);
                    [m,n]=size(no_zero);
                    if m~=0
                        B=D(i-2:i+2,j-2:j+2);
                    else
                         no_zero=find(D(i-3:i+3,j-3:j+3)~=0);
                         B=D(i-3:i+3,j-3:j+3);
                    end%以上是对没有被填充的点进行处理
                end
                cut_picture2(i-3,j-3)=median(B(no_zero));%原没有被插值成功的像素点，以周围非零像素的中位数代替；
            else
               cut_picture2(i-3,j-3)=cut_picture(i-3,j-3);
            end
        end
    end
    OutPicture(:,:,k)=cut_picture2;
end%这里是对图像三通道进行处理的，便于得到彩色图像

fig = figure(1);
imshow(OutPicture,'border','tight','initialmagnification','fit');
set(fig,'Position',[100,100,w,h]);
axis normal;
imwrite(OutPicture,[dir_fig num2str(t,'%.1f') '.jpg']);

end