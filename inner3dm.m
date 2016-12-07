function y = inner3dm(x)

[S1,S2,S3]    = size(x);
[S4,S5,S6] = size(x{1});
y       = zeros(S1,S2,S3,S4,S5,S6);

for i = 1:S1
    for j = 1:S2
        for k = 1:S3
            tmp          = x{i,j};
            try   tmp = tmp(1:S4,1:S5,1:S6);
            catch tmp = bulk(tmp,S4,S5,S6);
                  tmp = tmp(1:S4,1:S5,1:S6);
            end
            y(i,j,k,:,:,:) = tmp;
        end
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
