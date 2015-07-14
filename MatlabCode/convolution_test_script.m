% Script to test associativity of convolution in 3 dimensions
A = zeros(5,5,7);
B = zeros(3,3);
C = zeros(4,1);
test1 = zeros(5,5,7);
test2 = zeros(5,5,7);
test3 = zeros(5,5,7);

for i = 1:length(A)
    for j = 1:5
        for k = 1:5
            A(j,k,i) = random('unif',1,10);
        end
    end
end

for i = 1:length(B)
    for j = 1:length(B)
        B(i,j) = random('unif',1,10);
    end
end

for i = 1:length(C)
    C(i) = random('unif',1,10);
end

% test sequence 1
for i = 1:length(A)
    test1(:,:,i) = conv2(A(:,:,i),B,'same');
end
for i = 1:5
    for j = 1:5
        test1(i,j,:) = conv(squeeze(test1(i,j,:)),C,'same');
    end
end
disp(test1(:,:,2));

% test sequence 2
for i = 1:5
    for j = 1:5
        test2(i,j,:) = conv(squeeze(A(i,j,:)),C,'same');
    end
end
for i = 1:length(A)
    test2(:,:,i) = conv2(test2(:,:,i),B,'same');
end
disp(test2(:,:,2));