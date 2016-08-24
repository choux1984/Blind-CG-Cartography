function [Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset
%%
%This function returns the temperature data per month in H for every city
%the mean of the city in M and the altitude in Alt
%The subscript n and o denotes the old and new measurments
pth=which('readTemperatureDataset');
folder=[pth(1:end-(length('readTemperatureDataset')+2))  'TemperatureDataset/'];
file=strcat(folder,'temp_swiss1960');
Do = readtable(file,'Delimiter','\t','ReadVariableNames',false);
Ho=table2array(Do(:,5:end-1));
Mo=table2array(Do(:,end));
Alto=table2array(Do(:,2));
file=strcat(folder,'temp_swiss');
Dn = readtable(file,'Delimiter','\t','ReadVariableNames',false);
Hn=table2array(Dn(:,5:end-1));
Mn=table2array(Dn(:,end));
Altn=table2array(Dn(:,2));
end

