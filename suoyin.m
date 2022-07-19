% function [idxD,idxM]=suoyin(D,M,config,num)
% %num表示分割块数因子
% DD=double(D');
% MM=double(M');
% 
% DD_cell=cell(num,num);
% MM_cell=cell(num,num);
% shang1=zeros(num,num);
% shang2=zeros(num,num);
% index_DD=cell(num,num);
% index_MM=cell(num,num);
% 
% border1=[max(DD(:,1)),max(DD(:,2)),max(DD(:,3))];
% border11=[min(DD(:,1)),min(DD(:,2)),min(DD(:,3))];
% border2=[max(MM(:,1)),max(MM(:,2)),max(MM(:,3))];
% border22=[min(MM(:,1)),min(MM(:,2)),min(MM(:,3))];
% step1=(max((border1(1)-border11(1)),(border2(1)-border22(1))));
% step2=(max((border1(2)-border11(2)),(border2(2)-border22(2))));
% step3=(max((border1(3)-border11(3)),(border2(3)-border22(3))));
% % border_max=[max(border1(1),border2(1)),max(border1(2),border2(2)),max(border1(3),border2(3))];
% % border_min=[min(border11(1),border22(1)),min(border11(2),border22(2)),min(border11(3),border22(3))];
% size1=[step1/num,step2/num,step3/num];
% % step=(border1(1)-border11(1))/6;
% for i=1:num
%     for j=1:num
%         [index1,~,~]=find(((border11(1)+size1(1)*(i-1))<=DD(:,1)&DD(:,1)<=(border11(1)+size1(1)*i)) ...
%                      & ((border11(2)+size1(2)*(j-1))<=DD(:,2)&DD(:,2)<=(border11(2)+size1(2)*j)));
%         DD_cell{i,j}=DD(index1,:);
%         index_DD{i,j}=index1;
%         if(DD_cell{i,j}~=0)
%             shang1(i,j)=entropy(DD_cell{i,j});
%         else
%             shang1(i,j)=2;
%         end
%         [index2,~,~]=find(((border22(1)+size1(1)*(i-1))<=MM(:,1)&MM(:,1)<=(border22(1)+size1(1)*i)) ...
%                      & ((border22(2)+size1(2)*(j-1))<=MM(:,2)&MM(:,2)<=(border22(2)+size1(2)*j)));
%         MM_cell{i,j}=MM(index2,:);
%         index_MM{i,j}=index2;
%         if(MM_cell{i,j}~=0)
%             shang2(i,j)=entropy(MM_cell{i,j});
%         else
%             shang2(i,j)=2;
%         end
%     end
% end
% flag=1;
% % DD_cell
% % MM_cell
% % shang1
% % shang2
% shang1_arr=reshape(shang1,1,num*num);
% shang2_arr=reshape(shang2,1,num*num);
% for po=1:(num*num-1)
%     if abs(shang1_arr(po)-shang1_arr(po+1))<=0.00001
%         shang1_arr(po+1)=2;
%     end
%     if abs(shang2_arr(po)-shang2_arr(po+1))<=0.00001
%         shang2_arr(po+1)=2;
%     end
% end
% shang1=reshape(shang1_arr,num,num);
% shang2=reshape(shang2_arr,num,num);
% for i=1:num
%     for j=1:num
%         for ii=1:num
%             for jj=1:num
%                 if((shang1(i,j)~=2)&&(shang2(ii,jj)~=2)&&(length(cell2mat(DD_cell(i,j)))>=64)&&(length(cell2mat(MM_cell(ii,jj)))>=64))
%                     en{flag}=abs(shang1(i,j)-shang2(ii,jj));
%                     en_flag{flag}={i,j,ii,jj};
%                     flag=flag+1;
%                 end
%             end
%         end
%     end
% end
% [~,order]=sort(cell2mat(en));
% for i=1:length(order)
%     new_flag{i}=en_flag{order(i)};
% end
% idxD=[];
% idxM=[];
% for i=1:4
%     pp=[64,64,0,0];
%     uu=[cell2mat(new_flag{i}(1)),cell2mat(new_flag{i}(2)),cell2mat(new_flag{i}(3)),cell2mat(new_flag{i}(4))];
%     DD_part=DD_cell{uu(1),uu(2)};
%     MM_part=MM_cell{uu(3),uu(4)};
%     index1_new=index_DD{uu(1),uu(2)};
%     index2_new=index_MM{uu(3),uu(4)};
%     idx1 = randsample(size(DD_part,1), pp(i));
%     idx2 = randsample(size(MM_part,1), pp(i));
%     idxD=[idxD;index1_new(idx1)];
%     idxM=[idxM;index2_new(idx2)];
% end



function [idxD,idxM]=suoyin(D,M,config,num)
% %num表示分割块数因子
mark=64;%表示采样点数
DD=double(D');
MM=double(M');
enen=-1000*ones(num^2,num^2);
DD_cell=cell(num,num);
MM_cell=cell(num,num);
shang1=zeros(num,num);
shang2=zeros(num,num);
index_DD=cell(num,num);
index_MM=cell(num,num);

border1=[max(DD(:,1)),max(DD(:,2)),max(DD(:,3))];
border11=[min(DD(:,1)),min(DD(:,2)),min(DD(:,3))];
border2=[max(MM(:,1)),max(MM(:,2)),max(MM(:,3))];
border22=[min(MM(:,1)),min(MM(:,2)),min(MM(:,3))];
% step11=border1(1)-border11(1);
% step12=border1(2)-border11(2);
% step13=border1(3)-border11(3);
% step21=border2(1)-border22(1);
% step22=border2(2)-border22(2);
% step23=border2(3)-border22(3);
% step1=[step11,step12,step13];
% step2=[step21,step22,step23];
step1=border1-border11;
step2=border2-border22;
[~,order1]=sort(step1);
[~,order2]=sort(step2);
% border_max=[max(border1(1),border2(1)),max(border1(2),border2(2)),max(border1(3),border2(3))];
% border_min=[min(border11(1),border22(1)),min(border11(2),border22(2)),min(border11(3),border22(3))];
size1=step1(order1)/num;
size2=step2(order2)/num;
border_D=border11(order1);
border_M=border22(order2);
win=max(size1(3),size2(3))*num*1.2;
% step=(border1(1)-border11(1))/6;
for i=1:num
    for j=1:num
%         flll=(border_D(3)+size1(3)*(i-1)+win);
%         fll2=(border_D(3)+size1(3)*(i-1));
        [index1,~,~]=find(((border_D(3)+size1(3)*(i-1))<=DD(:,3)&DD(:,3)<=(border_D(3)+size1(3)*(i-1)+win)) ...
                     & ((border_D(2)+size1(2)*(j-1))<=DD(:,2)&DD(:,2)<=((border_D(2)+size1(2)*(j-1))+win)));
        DD_cell{i,j}=DD(index1,:);
        index_DD{i,j}=index1;
        if(DD_cell{i,j}~=0)
            shang1(i,j)=entropy(DD_cell{i,j});
        else
            shang1(i,j)=-1000;
        end
        [index2,~,~]=find(((border_M(3)+size2(3)*(i-1))<=MM(:,3)&MM(:,3)<=(border_M(3)+size2(3)*(i-1)+win)) ...
                     & ((border_M(2)+size2(2)*(j-1))<=MM(:,2)&MM(:,2)<=((border_M(2)+size2(2)*(j-1))+win)));
        MM_cell{i,j}=MM(index2,:);
        index_MM{i,j}=index2;
        if(MM_cell{i,j}~=0)
            shang2(i,j)=entropy(MM_cell{i,j});
        else
            shang2(i,j)=-1000;
        end
    end
end
flag=0;
% DD_cell
% MM_cell
% shang1
% shang2
shang1_arr=reshape(shang1,1,num*num);
shang2_arr=reshape(shang2,1,num*num);
shang1_diff=shang1_arr;
shang2_diff=shang2_arr;
for po=1:(num*num-1)
    if abs(shang1_arr(po)-shang1_arr(po+1))<=0.0001
        shang1_diff(po+1)=-1000;
    end
    if abs(shang2_arr(po)-shang2_arr(po+1))<=0.0001
        shang2_diff(po+1)=-1000;
    end
end
for po=1:(num*num-num)
    if abs(shang1_arr(po)-shang1_arr(po+num))<=0.0001
        shang1_diff(po+num)=-1000;
    end
    if abs(shang2_arr(po)-shang2_arr(po+num))<=0.0001
        shang2_diff(po+num)=-1000;
    end
end
shang1=reshape(shang1_diff,num,num);
shang2=reshape(shang2_diff,num,num);
for i=1:num
    for j=1:num
        for ii=1:num
            for jj=1:num
                if(shang1(i,j)~=-1000)&&(shang2(ii,jj)~=-1000)
                    if(length(cell2mat(DD_cell(i,j)))>=mark*1.5)
                        if(length(cell2mat(MM_cell(ii,jj)))>=mark*1.5)
                            enen(j+num*(i-1),jj+num*(ii-1))=exp(-abs(shang1(i,j)-shang2(ii,jj)));
                            flag=flag+1;
%                         else; flag=0;
                        end
%                     else; flag=0;
                    end
%                 else; flag=0;
                end
                enen_flag{j+num*(i-1),jj+num*(ii-1)}={i,j,ii,jj};
            end
        end
    end
end
map=solveSDP3(enen,config);
[w,d]=find(map==1);
idxD=[];
idxM=[];
pp=[mark/2,mark/2];
if length(w)<length(pp)
    fprintf('分块不合理,采用随机采样 \n');
    idxD = randsample(size(D,2), mark);
    idxM = randsample(size(M,2), mark);
else
    for i=1:length(pp)
        new_flag{i}= enen_flag{w(i),d(i)};
        uu=[cell2mat(new_flag{i}(1)),cell2mat(new_flag{i}(2)),cell2mat(new_flag{i}(3)),cell2mat(new_flag{i}(4))];
        DD_part=DD_cell{uu(1),uu(2)};
        MM_part=MM_cell{uu(3),uu(4)};
        index1_new=index_DD{uu(1),uu(2)};
        index2_new=index_MM{uu(3),uu(4)};
        idx1 = randsample(size(DD_part,1), pp(i));
        idx2 = randsample(size(MM_part,1), pp(i));
        idxD=[idxD;index1_new(idx1)];
        idxM=[idxM;index2_new(idx2)];
    end
end







% function [idxD,idxM]=suoyin(D,M,num)
% %num表示分割块数因子
% mark=64;%表示采样点数
% DD=double(D');
% MM=double(M');
% idxD=[];
% idxM=[];
% enen=nan*ones(num^2,num^2);
% en={};
% DD_cell=cell(num,num);
% MM_cell=cell(num,num);
% shang1=zeros(num,num);
% shang2=zeros(num,num);
% index_DD=cell(num,num);
% index_MM=cell(num,num);
% 
% border1=[max(DD(:,1)),max(DD(:,2)),max(DD(:,3))];
% border11=[min(DD(:,1)),min(DD(:,2)),min(DD(:,3))];
% border2=[max(MM(:,1)),max(MM(:,2)),max(MM(:,3))];
% border22=[min(MM(:,1)),min(MM(:,2)),min(MM(:,3))];
% step1=(max((border1(1)-border11(1)),(border2(1)-border22(1))));
% step2=(max((border1(2)-border11(2)),(border2(2)-border22(2))));
% step3=(max((border1(3)-border11(3)),(border2(3)-border22(3))));
% % border_max=[max(border1(1),border2(1)),max(border1(2),border2(2)),max(border1(3),border2(3))];
% % border_min=[min(border11(1),border22(1)),min(border11(2),border22(2)),min(border11(3),border22(3))];
% size1=[step1/num,step2/num,step3/num];
% % step=(border1(1)-border11(1))/6;
% for i=1:num
%     for j=1:num
%         [index1,~,~]=find(((border11(1)+size1(1)*(i-1))<=DD(:,1)&DD(:,1)<=(border11(1)+size1(1)*i)) ...
%                      & ((border11(2)+size1(2)*(j-1))<=DD(:,2)&DD(:,2)<=(border11(2)+size1(2)*j)));
%         DD_cell{i,j}=DD(index1,:);
%         index_DD{i,j}=index1;
%         if(DD_cell{i,j}~=0)
%             shang1(i,j)=entropy(DD_cell{i,j});
%         else
%             shang1(i,j)=nan;
%         end
%         [index2,~,~]=find(((border22(1)+size1(1)*(i-1))<=MM(:,1)&MM(:,1)<=(border22(1)+size1(1)*i)) ...
%                      & ((border22(2)+size1(2)*(j-1))<=MM(:,2)&MM(:,2)<=(border22(2)+size1(2)*j)));
%         MM_cell{i,j}=MM(index2,:);
%         index_MM{i,j}=index2;
%         if(MM_cell{i,j}~=0)
%             shang2(i,j)=entropy(MM_cell{i,j});
%         else
%             shang2(i,j)=nan;
%         end
%     end
% end
% flag=1;
% % DD_cell
% % MM_cell
% % shang1
% % shang2
% shang1_arr=reshape(shang1,1,num*num);
% shang2_arr=reshape(shang2,1,num*num);
% shang1_diff=shang1_arr;
% shang2_diff=shang2_arr;
% for po=1:(num*num-1)
%     if abs(shang1_arr(po)-shang1_arr(po+1))<=0.1
%         shang1_diff(po+1)=nan;
%     end
%     if abs(shang2_arr(po)-shang2_arr(po+1))<=0.1
%         shang2_diff(po+1)=nan;
%     end
% end
% for po=1:(num*num-num)
%     if abs(shang1_arr(po)-shang1_arr(po+num))<=0.1
%         shang1_diff(po+num)=nan;
%     end
%     if abs(shang2_arr(po)-shang2_arr(po+num))<=0.1
%         shang2_diff(po+num)=nan;
%     end
% end
% shang1=reshape(shang1_diff,num,num);
% shang2=reshape(shang2_diff,num,num);
% for i=1:num
%     for j=1:num
%         for ii=1:num
%             for jj=1:num
%                 if(~isnan(shang1(i,j)))&&(~isnan(shang2(ii,jj)))&&(length(cell2mat(DD_cell(i,j)))>=max(length(D)/100,mark))&&(length(cell2mat(MM_cell(ii,jj)))>=max(length(D)/100,mark))
%                     enen(j+num*(i-1),jj+num*(ii-1))=abs(shang1(i,j)-shang2(ii,jj));
%                     en{flag}=abs(shang1(i,j)-shang2(ii,jj));
%                     en_flag{flag}={i,j,ii,jj};
%                     flag=flag+1;
%                 end
%             end
%         end
%     end
% end
% % enen
% [~,order]=sort(cell2mat(en));
% if isempty(en)
%     fprintf('初始点云集太小,采用随机采样 \n');
%     idxD = randsample(size(D,2), mark);
%     idxM = randsample(size(M,2), mark);
% else
%     for i=1:length(order)
%         new_flag{i}=en_flag{order(i)};
%     end
% 
%     for i=1:4
%         pp=[mark/2,mark/2,0,0];
%         uu=[cell2mat(new_flag{i}(1)),cell2mat(new_flag{i}(2)),cell2mat(new_flag{i}(3)),cell2mat(new_flag{i}(4))];
%         DD_part=DD_cell{uu(1),uu(2)};
%         MM_part=MM_cell{uu(3),uu(4)};
%         index1_new=index_DD{uu(1),uu(2)};
%         index2_new=index_MM{uu(3),uu(4)};
%         idx1 = unique(randsample(size(DD_part,1), pp(i)));
%         idx2 = unique(randsample(size(MM_part,1), pp(i)));
%         idxD=[idxD;index1_new(idx1)];
%         idxM=[idxM;index2_new(idx2)];
%     end
% end
% % xx=entropy(MM_cell{2,2})