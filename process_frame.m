function [] = process_frame(dot,t,i_frame)
global v
%G0: 0,0; 600,0; 600,0; 600,600
%G1: 446,434; 820,434; 621,617; 1151,617
%G2: 

dir_fig = 'C:\Users\HLiang\Desktop\Case study LP\VideoFrames\';

f = read(v,i_frame);
% dot=[446,434; 621,617; 1151,617;820,434];       %ȡ�ĸ��㣬���������ϣ����£�����,����
y=[dot(1,1),dot(2,1),dot(3,1),dot(4,1)];        %�ĸ�ԭ����
x=[dot(1,2),dot(2,2),dot(3,2),dot(4,2)];
h = 600; w = 600;
%�������µĶ��㣬��ȡ�ľ���,Ҳ����������������״
%�����ԭͼ���Ǿ��Σ���ͼ���Ǵ�dot��ȡ�õĵ���ɵ������ı���.:)
Y=[1,1,w,w];     
X=[1,h,h,1];

B=[X(1),Y(1),X(2),Y(2),X(3),Y(3),X(4),Y(4)]';   %�任����ĸ����㣬�����ұߵ�ֵ
%�����ⷽ���飬���̵�ϵ��
A=[x(1),y(1),1,0,0,0,-X(1)*x(1),-X(1)*y(1);             
   0,0,0,x(1),y(1),1,-Y(1)*x(1),-Y(1)*y(1);
   x(2),y(2),1,0,0,0,-X(2)*x(2),-X(2)*y(2);
   0,0,0,x(2),y(2),1,-Y(2)*x(2),-Y(2)*y(2);
   x(3),y(3),1,0,0,0,-X(3)*x(3),-X(3)*y(3);
   0,0 ,0,x(3),y(3),1,-Y(3)*x(3),-Y(3)*y(3);
   x(4),y(4),1,0,0,0,-X(4)*x(4),-X(4)*y(4);
   0,0,0,x(4),y(4),1,-Y(4)*x(4),-Y(4)*y(4)];%���任���������ʽ
cut_picture=zeros(h,w);
OutPicture=zeros(h,w,3);
OutPicture=cast(OutPicture,'uint8');%���ͱ任Ϊuint8
fa=inv(A)*B;
fa=reshape([fa;1],[3,3]);%fa���Ǳ任����
cut=zeros(h+1,w+1);
BW=roipoly(f,y,x);%������ѡ������Ķ�ֵ��ͼ��cut
[coordinate_x,coordinate_y]=find(BW==1);%���ѡ�����������λ��
coordinate=[coordinate_x,coordinate_y,ones(length(coordinate_y),1)];%����z=1ֵ�Ա������߱任
cut=coordinate*fa;%��������仯
cut(:,1)=cut(:,1)./cut(:,3);
cut(:,2)=cut(:,2)./cut(:,3);%���ͼ���϶�Ӧ��x��yֵ
cut=floor(cut+eps);
cut(find(cut(:,1:2)==0))=1;%�߽�Ϊ��ĵ��Ϊ1��
cut_coordinate=cut(:,1:2);%��ñ任�������λ��
for k=1:3
    img=f(:,:,k);
    for i=1:length(coordinate)
        cut_picture(cut_coordinate(i,1),cut_coordinate(i,2))=img(coordinate(i,1),coordinate(i,2));
    end%����Ӧ��������ض�Ӧ
    cut_picture=cast(cut_picture,'uint8');
    B=zeros(5,5);
    D = padarray(cut_picture,[3 3],'replicate','both');%�������
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
                    end%�����Ƕ�û�б����ĵ���д���
                end
                cut_picture2(i-3,j-3)=median(B(no_zero));%ԭû�б���ֵ�ɹ������ص㣬����Χ�������ص���λ�����棻
            else
               cut_picture2(i-3,j-3)=cut_picture(i-3,j-3);
            end
        end
    end
    OutPicture(:,:,k)=cut_picture2;
end%�����Ƕ�ͼ����ͨ�����д���ģ����ڵõ���ɫͼ��

fig = figure(1);
imshow(OutPicture,'border','tight','initialmagnification','fit');
set(fig,'Position',[100,100,w,h]);
axis normal;
imwrite(OutPicture,[dir_fig num2str(t,'%.1f') '.jpg']);

end