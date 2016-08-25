function [C] = extractData()

% fileID = fopen('C:\Users\ioann006\Desktop\amazon-meta.txt');
% formatSpec='%s %f %f %f %s';
% N = 5;
% k = 0;
% while ~feof(fileID)
% 	k = k+1;
% 	C = textscan(fileID,formatSpec,N,'CommentStyle','\n\n','Delimiter','\t');
% 	figure(k)
% 	scatter(C{2},C{3})
% 	title(['Temperature and Humidity, Block ',num2str(k)])
% end


% fid = fopen('C:\Users\ioann006\Desktop\amazon-meta.txt', 'rt');
% formatString = repmat('%f,', 1, 1200);
% formatString (end) = [];
% allData = zeros(0, 1200);
% while ~feof(fid)
%   data = textscan(fid, formatString, 100000);
%   allData = [allData; [data{:}]];
% end
% fclose(fid);


fid = fopen('C:\Users\ioann006\Desktop\amazon-meta.txt');
C = cell(560000,10);
% while ~feof(fid)
C  = textscan(fid,'%s','Delimiter','\n');
% S=char(Out{1});
% Sini=strsplit(Out{1},' ');
% if(strcmp(Sini(1),'Id:'))
%     
%     id=str2double(Sini(2));
%     id=id+1;
%     Out  = textscan(fid,'%s',1,'Delimiter','\n');
%     S=char(Out{1});
%     %keyboard
%     Sini=strsplit(S,' ');
%     if(strcmp(Sini(1),'ASIN:'))
%     C(id,1)=(Sini(2));
%     Out  = textscan(fid,'%s',1,'\n','Delimiter','\n');
%     
%     end
%     if id==100
%     break;
%     end
% end
% end
end


