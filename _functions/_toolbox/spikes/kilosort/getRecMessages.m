function [msgs,ts] = getRecMessages(fnm,fnt)

% function [msgs, ts] = getRecMsgTimes(fnm,fnt)
% 
% this function reads a text.npy file and timestamps.npy file
% generated by openEphys in flat binary format. displays the
% extracted events to the command window.
%   
% input: [fnm, fnt] filenames for text and timestamps as strings
% output: [msgs, ts] message text, timestamps in samples
%
% TO DO: check that the timestamps are offset from pressing play vs record!
    
% open the file
fid = fopen(fnm,'r');
msg = textscan(fid,'%s', 'delimiter',{'\n'});
fclose(fid);

if any(size(msg{1}) > 1)
    % split messages by '0' chars (these are null)
    msgs = strsplit(msg{1}{2},'\0')';

    % remove empties
    msgs = msgs(~cellfun('isempty',msgs));
else
    msgs = [];
end

% get timestamps
ts = double(readNPY(fnt));

fprintf('\n');
fprintf('\tTimestamp Message\n\t-----------------\n');
for i = 1:length(msgs)
    fprintf('\t%d %d %s\n',i,ts(i),msgs{i})
end
fprintf('\n');