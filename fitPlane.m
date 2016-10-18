function [n,b,bestn] = fitPlane(x,ns,nsample)
N = size(x,2);
thres = 5e-4;
beste = -1;
bestn = [];
for i = 1:nsample
    idx = randi(N,1,ns);
    px = x(:,idx);
    [n,b] = computePlane(px);
    if sum(abs(n))~=0
        err = abs(n.'*x+b);
        e = sum(err<thres);
        if e>beste
            beste = e;
            bestn = err<thres;
        end
    end
end
[n,b] = computePlane(x(:,bestn));
err = abs(n.'*x+b);
bestn = err<3e-3;
end

function [n,b] = computePlane(x)
n = [0 0 0].';
b = 0;
mx = mean(x,2);
cx = x - repmat(mx,1,size(x,2));
covx = cx*cx.';
[Q,l] = eig(covx);
[val,idx] = min(abs(diag(l)));
if val < 0.01
    n = Q(:,idx);
    b = -n.'*mx;
end
end