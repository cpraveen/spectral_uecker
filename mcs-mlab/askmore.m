% little function to ask user for stop or continue, with default 
function o=askmore(i); 
if(i==1) reply= input('More? Y/n: ', 's'); 
else reply= input('More? y/N: ', 's'); end
if isempty(reply) o=i; elseif(i==1) o=0; else o=1; end 
end