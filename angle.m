 function du = angle(M,p1,p2,p3,config)
 % side(1) = sqrt(power(a(1) - b(1),2)+pow(a(2) - b(2),2) + pow(a(3) - b(3),2)); 
 % side(2) = sqrt(power(a(1) - c(1),2)+pow(a(2) - c(2),2) + pow(a(3) - c(3),2));
 % side(3) = sqrt(power(c(1) - b(1),2)+pow(c(2) - b(2),2) + pow(c(3) - b(3),2)); 
 side(1) = abs( distance(M,p1,p2));
 side(2) = abs( distance(M,p2,p3));
 side(3) = abs( distance(M,p1,p3));
 side_side=[M(:,p1),M(:,p2);M(:,p2),M(:,p3);M(:,p1),M(:,p3)];
 %% 条件基于三角形的基本定律以及保证边长不过小
 if((side(1)+side(2)<=side(3)) || (side(1)+side(3)<=side(2)) || (side(2)+side(3)<=side(1)) ...
         && side(1) >= config.pairDistThresh && side(2) >= config.pairDistThresh && side(3) >= config.pairDistThresh) 
     du(1)=0;
     du(2)=0;
 else
     [~,order]=sort(side);
     side_flag=side(order);
     side_side_flag=side_side(order,:);
     du(1)=side_side_flag(1,1)'*side_side_flag(1,2)/(side_flag(1)+side_flag(2));
     du(2)=side_side_flag(2,1)'*side_side_flag(2,2)/(side_flag(2)+side_flag(3));
 end
 end