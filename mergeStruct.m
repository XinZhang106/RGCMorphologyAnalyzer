function new_struct = mergeStruct(stru1, stru2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
f1 = fieldnames(stru1);
f2 = fieldnames(stru2);
new_struct = stru1; % Initialize new_struct with the contents of stru1
for i = 1:length(f2)
    if (~ismember(f2{i}, f1))
        new_struct.(f2{i}) = stru2.(f2{i}); % Merge fields from stru2 into new_struct
    else
        new_struct.(f2{i}){end+1} = new_struct.(f2{i}); %append to the field with the same name
    end
end
end