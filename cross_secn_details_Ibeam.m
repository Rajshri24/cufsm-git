
function [span, d, bf, tf, tw] = cross_secn_details_Ibeam()

span = 480;

%Provide the W-section's name below and it takes the cross-sectional
%properties from the aisc shapes database
section = 'W24X68';

[num,text,raw] = xlsread('aisc-shapes-database-v15.0.xlsx', 2);
p = strcmp(section,text(:,3));
rowNum = find(p==1); 

d = num(rowNum-1, 3);
bf = num(rowNum-1, 8);
tf = num(rowNum-1, 16);
tw = num(rowNum-1, 13);

end