function  A = mfp1p2gen(n,prob1,prob2,iterations)
%generates k iterations of mfp1p2 model of size nxn with parameters prob1 and prob2

p1 = prob1;
p2 = prob2;
k = iterations;

sz = n;
sqsz = 2.^(1:k);
sqsz = repelem(sqsz,2);
sqsz(k+1:end) = [];

initmat = ones(sz);
probmat = ones(sz);

% probmat(1:end/2,:) = p1;
% probmat((end/2)+1:end,:) = p2;

rowvec = sz;
colvec = sz;

counter = 0;

for x = sqsz
    counter = counter + 1;

    if mod(counter,2) == 0
        pxsz = floor(sz./x); 

        numboxC = floor(sz / pxsz); % vertical cut
        colvec = pxsz * ones(1, numboxC);

        probcell = mat2cell(probmat,rowvec,colvec);

        for i = 1:length(probcell(:,1))
            for j = 1:2:length(probcell(1,:))-1
                left = cell2mat(probcell(i,j));
                right = cell2mat(probcell(i,j+1));

                half = rand;
                if half < 0.5
                    left = left.*p1;
                    right = right.*p2;
                else
                    left = left.*p2;
                    right = right.*p1;
                end

                probcell{i,j} = left;
                probcell{i,j+1} = right;

                
            end
        end

        probmat = cell2mat(probcell);

    else
        pxsz = floor(sz./x);

        numboxR = floor(sz / pxsz); % horizontal cut
        rowvec = pxsz * ones(1, numboxR);

        probcell = mat2cell(probmat,rowvec,colvec);

        for i = 1:2:length(probcell(:,1))-1
            for j = 1:length(probcell(1,:))
                top = cell2mat(probcell(i,j));
                bot = cell2mat(probcell(i+1,j));

                half = rand;
                if half < 0.5
                    top = top.*p1;
                    bot = bot.*p2;
                else
                    top = top.*p2;
                    bot = bot.*p1;
                end

                probcell{i,j} = top;
                probcell{i+1,j} = bot;
            end
        end

        probmat = cell2mat(probcell);

    end
end

initcell = mat2cell(initmat,rowvec,colvec);

for i = 1:length(probcell(:,1))
    for j = 1:length(probcell(1,:))
        len1 = size(probcell{i,j},1);
        len2 = size(probcell{i,j},2);
        temp = initcell{i,j};
        tempprob = probcell{i,j};
        num = ceil((len1*len2)*tempprob(1,1));
        mycount = 0;
        for v = 1:len1
            for w = 1:len2
                if mycount < num
                    temp(v,w) = 0;
                    mycount = mycount+1;
                end
            end
        end
        myvec = reshape(temp,1,[]);
        myvec = myvec(randperm(length(myvec)));
        temp = reshape(myvec,size(temp));
        initcell{i,j} = temp;
    end
end

A = cell2mat(initcell);