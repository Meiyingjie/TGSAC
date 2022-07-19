close all
dataset = 'synthetic';
method = 'SDRSAC';          % Without Correspondences


% Read configuration containing hyperparameters
config = readConfig(dataset); 

% Read Data
load(config.matPath);
a=pcread(config.plyPath);
b=pcread(config.plyPathB);
% figure;
% pcshow(a);
% hold on;
% pcshow(b);
bb=pcdownsample(b,'random',0.1);  
aa=pcdownsample(a,'random',0.1);  
D=bb.Location';
M=aa.Location';
DD=double(D');
MM=double(M');
num=10;

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
step1=(max((border1(1)-border11(1)),(border2(1)-border22(1))));
step2=(max((border1(2)-border11(2)),(border2(2)-border22(2))));
step3=(max((border1(3)-border11(3)),(border2(3)-border22(3))));
% border_max=[max(border1(1),border2(1)),max(border1(2),border2(2)),max(border1(3),border2(3))];
% border_min=[min(border11(1),border22(1)),min(border11(2),border22(2)),min(border11(3),border22(3))];
size1=[step1/num,step2/num,step3/num];
% step=(border1(1)-border11(1))/6;
for i=1:num
    for j=1:num
        [index1,~,~]=find(((border11(1)+size1(1)*(i-1))<=DD(:,1)&DD(:,1)<=(border11(1)+size1(1)*i)) ...
                     & ((border11(2)+size1(2)*(j-1))<=DD(:,2)&DD(:,2)<=(border11(2)+size1(2)*j)));
        DD_cell{i,j}=DD(index1,:);
        index_DD{i,j}=index1;
        if(DD_cell{i,j}~=0)
            shang1(i,j)=entropy(DD_cell{i,j});
        else
            shang1(i,j)=0;
        end
        [index2,~,~]=find(((border22(1)+size1(1)*(i-1))<=MM(:,1)&MM(:,1)<=(border22(1)+size1(1)*i)) ...
                     & ((border22(2)+size1(2)*(j-1))<=MM(:,2)&MM(:,2)<=(border22(2)+size1(2)*j)));
        MM_cell{i,j}=MM(index2,:);
        index_MM{i,j}=index2;
        if(MM_cell{i,j}~=0)
            shang2(i,j)=entropy(MM_cell{i,j});
        else
            shang2(i,j)=0;
        end
    end
end
flag=1;
for i=1:num
    for j=1:num
        for ii=1:num
            for jj=1:num
%                 length(MM_cell(i,j))
                if((shang1(i,j)~=0)&&(shang2(ii,jj)~=0)&&(length(cell2mat(DD_cell(i,j)))>=5)&&(length(cell2mat(MM_cell(ii,jj)))>=5))
                    en{flag}=abs(shang1(i,j)-shang2(ii,jj));
                    en_flag{flag}={i,j,ii,jj};
                    flag=flag+1;
                end
            end
        end
    end
end
[~,order]=sort(cell2mat(en));
for i=1:length(order)
    new_flag{i}=en_flag{order(i)};
end
idxD=[];
idxM=[];
for i=1:5
    pp=6-i;
    uu=[cell2mat(new_flag{i}(1)),cell2mat(new_flag{i}(2)),cell2mat(new_flag{i}(3)),cell2mat(new_flag{i}(4))];
    DD_part=DD_cell{uu(1),uu(2)};
    MM_part=MM_cell{uu(3),uu(4)};
    index1_new=index_DD{uu(1),uu(2)};
    index2_new=index_MM{uu(3),uu(4)};
    idx1 = randsample(size(DD_part,1), pp);
    idx2 = randsample(size(MM_part,1), pp);
    idxD=[idxD;index1_new(idx1)];
    idxM=[idxM;index2_new(idx2)];
end
% xx=entropy(MM_cell{2,2})