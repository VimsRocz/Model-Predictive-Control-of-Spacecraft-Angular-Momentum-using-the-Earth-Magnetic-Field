function S = skew(v)
%SKEW 3x1 -> 3x3 skew-symmetric matrix so that S*w = v x w
v = v(:);
S = [  0   -v(3)  v(2);
     v(3)   0   -v(1);
    -v(2)  v(1)   0  ];
end
