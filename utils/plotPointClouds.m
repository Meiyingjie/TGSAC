function plotPointClouds(M, B, markerm, markerb)
    if (nargin < 4)
        markerm = 'r.';
    end
    if (nargin < 3)
        markerb = 'b.';
    end
    %figure;
    if strcmp(markerb, '+k')
        plot3(M(1,:), M(2,:), M(3,:), markerm, B(1,:), B(2,:), B(3,:), markerb, 'MarkerSize',10);
    else
    plot3(M(1,:), M(2,:), M(3,:), markerm, B(1,:), B(2,:), B(3,:), markerb, 'MarkerSize',0.8);
    end
    axis off;    

end