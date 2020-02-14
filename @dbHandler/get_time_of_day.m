function time = get_time_of_day(~, key)
%GET_TIME_OF_DAY Takes a key and splits out total hours
%   That's it.
time = strsplit(key, '_');
time = time{3};
h = str2num(time(1:2));
m = str2num(time(3:4));
s = str2num(time(5:6));

time = h + m / 60 + s / 3600;
end

