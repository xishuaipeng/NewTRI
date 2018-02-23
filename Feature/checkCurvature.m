function checkCurvature(x,y,k)
      maxLength = length(x);
     % k = smooth(k,20);
      label = k;
      classes = 2;
      label(label<700)=1;
      label(label>700)=2;
      %label(label<1)=3;
      
      radio = jet( classes);
      figure; hold on;
      for i = 1: classes
          index = find(label == i);
          p_x = x(index);
          p_y = y(index);
          color = radio(i,:);
          scatter(p_x,p_y, 'MarkerFaceColor',color,'MarkerEdgeColor',color);  
      end

end