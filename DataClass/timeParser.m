function time = timeParser(date)
     % tranfer the time from 'MM:SS.FFF' or 'HH:MM:SS.FFF' to seconds
         colonIndex = find(date == ':');
         if length(colonIndex)==1
             dateVector = datevec(date, 'MM:SS.FFF');
         elseif length(colonIndex)==2
             dateVector = datevec(date, 'HH:MM:SS.FFF');
         else
                'Time can not formated to second'
                time = date;
                return 
         end
         time =  seconds(duration(dateVector(:,4),dateVector(:,5), dateVector(:,6)));     
end

