clear variables;

U = 17.99;

Nx = 10;
Ny = 10;
epsilon = 0.01;
Ec = 1;

R0 = 2;
dR = 2;

Ru = R0 * ones(Ny-1, Nx-1) + dR * rand(Ny-1, Nx-1);
Rv = R0 * ones(Ny, Nx) + dR * rand(Ny, Nx);

Iu = zeros(Ny-1, Nx-1);
Iv = zeros(Ny, Nx);
E = zeros(Ny-1, Nx);
nIu = length(Iu(:));
nIv = length(Iv(:));
nE = length(E(:));

nVar = nIu + nIv + nE;
figure('Color', 'k');

for t = 1 : 100
    A = zeros(nVar, nVar);
    b = zeros(nVar, 1);

    % loop for each joint, the current is divergence-free
    index = 1;
    for i = 1 : Ny-1
        up = i + length(Iu(:));
        down = up + 1;
        right = i;

        A(index, up) = 1;
        A(index, down) = -1;
        A(index, right) = -1;
        index = index + 1;

        for j = 2 : Nx-1
            up = i + (j-1) * Ny + length(Iu(:));
            down = up + 1;
            left = i + (j-2) * (Ny - 1);
            right = i + (j-1) * (Ny - 1);

            A(index, up) = 1;
            A(index, down) = -1;
            A(index, left) = 1;
            A(index, right) = -1;

            index = index + 1;
        end

        up = i + (Nx-1) * Ny + length(Iu(:));
        down = up + 1;
        left = i + (Nx-2) * (Ny - 1);

        A(index, up) = 1;
        A(index, down) = -1;
        A(index, left) = 1;
        index = index + 1;
    end

    % loop for each resistor, following Ohm's law
    % Iu = zeros(Ny-1, Nx-1);
    for i = 1 : Ny - 1
        for j = 1 : Nx - 1
            left = i + (j - 1) * (Ny - 1) + length(Iu(:)) + length(Iv(:));
            right = left + (Ny - 1);
            mid = i + (j - 1) * (Ny - 1);

            A(index, mid) = Ru(i, j);
            A(index, left) = -1;
            A(index, right) = 1;

            index = index + 1;
        end
    end

    % Iv = zeros(Ny, Nx);
    for j = 1 : Nx
        down = 1 + (j - 1) * (Ny - 1) + length(Iu(:)) + length(Iv(:));
        mid = 1 + (j - 1) * Ny + length(Iu(:));

        A(index, mid) = Rv(1, j);
        A(index, down) = 1;
        b(index) = U;
        index = index + 1;
    end

    for i = 2 : Ny-1
        for j = 1 : Nx
            up = i - 1 + (j - 1) * (Ny - 1) + length(Iu(:)) + length(Iv(:));
            down = up + 1;
            mid = i + (j - 1) * Ny + length(Iu(:));

            A(index, mid) = Rv(i, j);
            A(index, up) = -1;
            A(index, down) = 1;

            index = index + 1;
        end
    end

    for j = 1 : Nx
        up = Ny - 1 + (j - 1) * (Ny - 1) + length(Iu(:)) + length(Iv(:));
        mid = Ny + (j - 1) * Ny + length(Iu(:));  

        A(index, up) = -1;
        A(index, mid) = Rv(Ny, j);
        b(index) = 0;

        index = index + 1;
    end
    
    A = sparse(A);
    x = A \ b;

    Iu = abs(reshape(x(1:nIu), Ny-1, Nx-1));
    Iv = abs(reshape(x((nIu+1):(nIu+nIv)), Ny, Nx));
    E = reshape(x((nIu+nIv+1):end), Ny-1, Nx);

    minI = min([Iu(:);Iv(:)]);
    maxI = max([Iu(:);Iv(:)]);

    clf;
    for i = 1 : Ny - 1
        for j = 1 : Nx - 1
            y = Ny - i;
            x = j - 1;
            color = interp_color(Iu(i, j), minI, maxI);
            line([x, x+1], [y, y], 'Color', color, 'LineWidth', 2);
        end
    end

    for i = 1 : Ny
        for j = 1 : Nx
            x = j - 1;
            y = Ny - i;
            color = interp_color(Iv(i, j), minI, maxI);
            line([x, x], [y, y+1], 'Color', color, 'LineWidth', 2);
        end
    end

    axis equal tight;
    axis off;
    
    pause(0.01);
    

    dEu = Iu .* Ru;
    dEv = Iv .* Rv;

    breakIndexU = find(dEu >= Ec);
    breakIndexV = find(dEv >= Ec);
    
    if length(breakIndexU) + length(breakIndexV) == 0
        break;
    end
    
    for index = breakIndexU'
        i = mod(index, Ny-1);
        if eq(i, 0)
            i = Ny - 1;
        end
        j = 1 + (index - i) / (Ny - 1);
        Ru(i, j) = epsilon * Ru(i, j);
    end

    for index = breakIndexV'
        i = mod(index, Ny);
        if eq(i, 0)
            i = Ny;
        end
        j = 1 + (index - i) / Ny;
        Rv(i, j) = epsilon * Rv(i, j);
    end
end


