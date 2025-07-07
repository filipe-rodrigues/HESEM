function delta_Z = calcula_Delta_Z(k,z,W,medidas,de,Ysh,V,yykm,TAP,bbus,para)
    
    q_medidas = length(medidas(:,1));   % number of measurement
    tipo = medidas(:,2);                % Type
    delta_Z = zeros(q_medidas,1);
    
    if k == 2
        for i = 1:q_medidas
           if tipo(i) == 1
               delta_Z(i,1) = z(i)-1;
               
           elseif tipo(i) == 2
               delta_Z(i,1) = conj(z(i)) - Ysh(de(i));
               
           elseif tipo(i) == 3
               delta_Z(i,1) = conj(z(i)) - ( (1/TAP(de(i),para(i)))^2 * yykm(de(i),para(i)) - (1/TAP(de(i),para(i))) * yykm(de(i),para(i)) + 1j*bbus(de(i),para(i))) ;
               
           else
               delta_Z(i,1) = z(i) - ( (1/TAP(de(i),para(i)))^2 * yykm(de(i),para(i)) - (1/TAP(de(i),para(i))) * yykm(de(i),para(i)) + 1j*bbus(de(i),para(i))) ;
           end
        end
        
    else
        for i = 1:q_medidas
           if tipo(i) == 1
               delta_Z(i,1) = 0;
               
           elseif tipo(i) == 2
               delta_Z(i,1) = conj( z(i) ) * conj( W( de(i),k-1 ) ) - Ysh(de(i))*V( de(i),k-1 );
               
           elseif tipo(i) == 3
               delta_Z(i,1) = conj(z(i))*conj(W(de(i),k-1));
               
           else
               delta_Z(i,1) = 0;
           end
        end
        
    end
end

