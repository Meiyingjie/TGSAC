function [idx1,idx2]=TAB(p,q,num)
    xx1=tabulate(p);
    xx2=tabulate(q);
    [~,order1]=sort(xx1(:,3),'descend');
    [~,order2]=sort(xx2(:,3),'descend');
    roro1=xx1(:,1);
    roro2=xx2(:,1);
    roro11=roro1(order1);
    roro22=roro2(order2);
    idx1=roro11(1:num);
    idx2=roro22(1:num);
end