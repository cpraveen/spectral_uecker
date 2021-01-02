function o=ask(s,i); % ask user  
% setup for easy hitting return resp easy stopping 
reply= input(s, 's'); 
if isempty(reply) o=i; else o=reply; end 
end 