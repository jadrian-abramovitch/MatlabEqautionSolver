% laplace3d
%
% We are approximating the solution to the Laplace differential equation using a grid 
% of points to generate a system of equations to solve, in either 2 or 3 dimensions. We % will observe either insulated or Dirichlet boundary conditions and adjust our matrix % accordingly.  
%
%
% Parameters
% ==========
%    U The matrix on which we are solving. Specifies our boundary conditions.
%
% Return Values
% =============
%    U_out The matrix Uout is the solved matrix U.

function [U_out] = laplace3d( U )

% Error checking
% ==============
%
%Check that U is a matrix, none of the boundary points = -Inf, must be boundary values around entire array
if ndims(U)~=3
    throw( MException( 'MATLAB:invalid_argument', ...
        'the argument U is not a 3-dimensional matrix' ) );
end
%if (sum(U(1,:)) == -Inf) || (sum(U(end,:)) == -Inf) || (sum(U(:,1)) == -Inf) || (sum(U(:,end))== -Inf)
   % throw( MException( 'MATLAB:invalid_argument', ...
        %'the argument U does not have only boundary conditions at the edges' ) );    
%end

% Initialization
% ==============
%
% Get the size of the matrix, assign the output to the input matrix
[n_x, n_y, n_z] = size( U );
U_out = U;
% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
% Iterate through matrix U, assigning an increasing integer value to each -Inf in the u_to_w matrix
% Assign the co-ordinated of each -Inf to w_to_u
u_to_w = zeros( n_x, n_y, n_z );
w_to_u = zeros( 3, n_x * n_y* n_z );
m = 0;
for ix = 1:n_x
    for iy = 1:n_y
        for iz = 1:n_z
        if U(ix, iy, iz) == -Inf
            m = m + 1;
            u_to_w(ix, iy, iz) = m;
            w_to_u(:, m) = [ix, iy, iz]';
        end
        end
    end
end
% Creating and solving a system of linear equations
% =================================================
%
% Create M and b matricies, then solve Mw = b
% Need to check the values around point in question, check what kind of
% value it is and treat it appropriately
M = zeros( m, m);
b = zeros(m, 1); 
for i = 1:m
    c = w_to_u(:,i);
    for k = 1:6
        if k ==1
            p = c + [-1, 0, 0]';
        elseif k==2
            p = c + [1, 0, 0]';
        elseif k ==3
            p = c + [0, -1, 0]';
        elseif k==4
            p = c + [0, 1, 0]';
        elseif k ==5
            p = c + [0, 0, -1]';
        elseif k==6
            p = c + [0, 0, 1]';
        end
        if (U(p(1),p(2),p(3))) == -Inf
            M(i,i) = M(i,i) -1;
            j = u_to_w(p(1),p(2),p(3));
            M(i,j) = M(i,j) +1;
        elseif isnan(U(p(1),p(2),p(3))) == 1
        else
            M(i,i) = M(i,i)-1;
            b(i) = b(i)-U(p(1),p(2),p(3));
        end
    end
end
U = M\b;
% Substituting the values back into the matrix Uout
% =================================================
%
% Use the w_to_u matrix to convert back into Uout, then return the value to the user.

for i = 1:m
    point = w_to_u(:,i);
    U_out(point(1),point(2),point(3)) = U(i);
end    
end
