%A MATLAB script that algorithmically converts a matrix into row-echelon form without pivoting, then solves the systems of equations
function f = naive_gauss(A,b)

    %Check for square matrix
    [n, m] = size(A);
    if n == m

        %Check if size of A and b match
        [x, ~] = size(b);
        if n == x

            %Check whether A is a singular matrix using det function
            try
                det(A);
            catch
                error("Matrix [A] is a singular matrix")
            end

            %Combine matrix A and vector b
            C = horzcat(A,b);

            %Perform Naive Gaussian Elimination

            %Loop for each column to be eliminated
            for s = 1:n-1

                %Loop for each row
                for i = s:n-1 
                    factor = C(i+1, s) / C(s, s);

                    %Loop for each cell in the row
                    for j = s:n+1
                        C(i+1, j) = C(i+1, j) - (factor * C(s, j));
                    end
                end
            end

            %Initialize x
            x = zeros(n,1);

            %Calculate x
            for i = n:-1:1
                intermediateVal = C(i, n + 1);
                for j = n:-1:i
                    if i == n
                        intermediateVal = C(i, n + 1);
                    else                        
                        intermediateVal = intermediateVal - (C(i, j) * x(j));
                    end
                end
                %disp(intermediateVal);
                x(i) = intermediateVal / C(i, i);
                %disp(x);
            end
            f = x;
        else
            error("Matrix dimension mismatch")
        end
    else
        error("Matrix [A] is not a square matrix")
    end
end
