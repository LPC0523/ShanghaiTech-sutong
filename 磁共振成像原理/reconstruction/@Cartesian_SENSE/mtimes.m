function [res] = mtimes(a,b)

% a = encoding operator, 
% b = full_k_dUa : (Ny,Nx,Nc) -> (Ny,Nx) when applying EH
% b = image_space : (Ny,Nx) -> (Ny,Nx,Nc) when applying E

if a.adjoint % EH operation
    res_coils = zeros(size(b,1),size(b,2),size(a.C,3)); % init
    
    coils = a.C;
    for coil = 1:size(coils,3) % number of coils
        % Sampling -> FFT -> Coil
        aux_b = b(:,:,coil).* a.U;
        res_coils(:,:,coil) = ktoi(aux_b).*conj(coils(:,:,coil)); 
    end

    res = sum(res_coils,3);
    
    % Applying intensity correction
    res = res ./ a.coil_rss;
    res(isnan(res)) = 0;
    res(isinf(res)) = 0;

else % E operation
    res = zeros(size(b,1),size(b,2),size(a.C,3)); % init
    
    % Applying intensity correction
    b = b ./ a.coil_rss;
    b(isnan(b)) = 0; 
    b(isinf(b)) = 0;
    
    coils = a.C;
    for coil = 1:size(coils,3) % number of coils
        % Coil -> FFT -> Sampling
        res(:,:,coil) = itok(b.*coils(:,:,coil)).*a.U;
    end

end

