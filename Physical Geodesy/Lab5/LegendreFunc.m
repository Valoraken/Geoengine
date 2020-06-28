function [P] = LegendreFunc(l, t)

n = length(t);
P = zeros(l+1, l+1, n);
P(1,1,:) = ones(1, n);
P(2,1,:) = t;
P(2,2,:) = sqrt((1-t.^2));

for i = 2:1:l
    for j = 0:1:l
        if i == j
            P(i+1,j+1,:) = (2*i-1)*P(2,2,:).*P(i,j,:);
        else
            P(i+1,j+1,:) = ((2*i-1)*P(2,1,:).*P(i,j+1,:)-(i-1+j)*P(i-1,j+1,:))/(i-j);
        end
    end
end

end
