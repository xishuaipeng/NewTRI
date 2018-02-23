function kappa = Curvature(x,y)
    x_ = diff(x,1);
    x_ = x_(2:end);
    x__ = diff(x,2);
    X = x(3:end);
    
    
    y_ = diff(y,1);
    y_ = y_(2:end);
    y__ = diff(y,2);
    Y = y(3:end);
    
    kappa = abs(x_.* y__ - x__ .* y_)./(x_.^2 + y_.^2).^(3/2);
    kappa = [kappa(1);kappa(1);kappa];
    
    
     %num = length(x);
%      stride = 0;
%      width = 1;
%      height = 1;
%      minx = min(x(:));
%      maxx = max(x(:));
%      miny =  min(y(:));
%      maxy =  max(y(:));
%      max_range = max((maxx -minx ),(maxy-miny));
%      x = stride+(width-2*stride) *(x-minx)/max_range;%(maxx-minx);
%      y = stride+(height-2*stride) *(y-miny)/max_range;%(maxy-miny);
%      p = polyfit(x , y, 2);
%      y_= polyval(p,x);
%      %compute the curvature
%      pp1 = 2 * p(1) * y_ + p(2);
%      pp2 = 2 * p(1);
%      kappa = abs(pp2) ./ sqrt((1+pp1.^2).^3);
end

