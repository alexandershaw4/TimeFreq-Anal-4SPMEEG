function y = inner3d(x)

[S1,S2]    = size(x);
[S3,S4,S5] = size(x{1});
y       = zeros(S1,S2,S3,S4,S5);

for i = 1:S1
    for j = 1:S2
        tmp          = x{i,j};
        try   tmp = tmp(1:S3,1:S4,1:S5);
        catch tmp = bulk(tmp,S3,S4,S5);
              tmp = tmp(1:S3,1:S4,1:S5);  
        end
        y(i,j,:,:,:) = tmp;
    end
end

end

function o = bulk(tmp,s1,s2,s3)

i = find(size(tmp)~=[s1 s2 s3]);
s = zeros(s1, s2, s3);
d = abs(size(tmp,i)-size(s,i));

switch i
    case 3
        tmp(:,:,end+1:end+d)=zeros(s1,s2,d);
end



o = tmp;
end
