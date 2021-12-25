%Written & Copyright: V. Kumar, vkumar@utep.edu
data = struct('fn',fn); 


%-----------------------------------------------
api='https://classes.mu2com.com/matlab.php'; 
url = sprintf('%s?uid=%s&pw=%s&fn=%s',api,uid,pin,fn); 
filetext = fileread([fn '.m']); data.filetext = filetext; 
msg = webwrite(url,'data',matlab.internal.webservices.toJSON(data));
%msg = webwrite(url,'data',jsonencode(data));
fprintf('-------------\n%s\n-------------\n',strtrim(msg));

