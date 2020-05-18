function [TIME] = get_time_of_day(~, key)
%GET_TIME_OF_DAY Takes a key and splits out total hours
%   That's it.
key_split = strsplit(key, '_');
day = key_split{end-2};
time = key_split{end-1};
hour = time(1:2);
h = str2num(hour);
minute = time(3:4);
m = str2num(minute);
second = time(5:6);
s = str2num(second);

hours_tot = h + m / 60 + s / 3600;

year  = num2str(day(1:2));
month = num2str(day(3:4));
date  = num2str(day(5:6));
TIME = datetime([year '-' month '-' date ' ' hour ':' minute ':' second] ,'InputFormat','yy-MM-dd HH:mm:ss');
end

