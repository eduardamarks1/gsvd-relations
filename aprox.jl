using Printf
using LinearAlgebra


function Make_D1_D2(bl,br,wl,wr,r, D1_til, D2_til)
    D1= [I(bl)        zeros(bl,r) zeros(bl,br);
         zeros(r,bl)  D1_til      zeros(r,br);
         zeros(wl,bl) zeros(wl,r) zeros(wl,br)]

    #D2 in our format
    D2= [zeros(wr,bl) zeros(wr,r) zeros(wr,br);
         zeros(r,bl)  D2_til      zeros(r,br);
         zeros(br,bl) zeros(br,r) I(br)]

    return D1, D2
end

function Make_A_B(bl,br,wl,wr,r,D1_til,D2_til,H,U,V)
    D1,D2 = Make_D1_D2(bl,br,wl,wr,r, D1_til, D2_til)

    A = U * D1 * H
    B = V * D2 * H

    return A, B
end

function Our_SVD(A,B)
    bl,br,wl,wr,r = wire_size(A,B)
    U, V, Q, D1, D2_julia, R0 = svd(A, B)

    P, D2_our = permutation(wr,r,br, D2_julia)
    V_til = V * P #Changing variable to fix line changes in D2

    #Selects D1_tilde and D2_tilde
    D1_til = D1[bl+1:bl+r, bl+1:bl+r]
    D2_til = D2_our[wr+1:wr+r, bl+1:bl+r]

    H = R0 * Q'

    return U, V_til, H, D1_til, D2_til, D1, D2_our, bl, br, wl, wr, r, P, V
end

function wire_size(A,B)
    m,.. = size(A)
    p,.. = size(B)
    .., .., .., D1, D2, .. = svd(A, B)
    ..,k_mais_l = size(D1)

    #Counting k: D1 is a m-by-(k_mais_l) diagonal matrix with 1s in the first k entries,
    k = 0
    for i in 1:min(m, k_mais_l)
        if D1[i, i] == 1
            k += 1
        else
            k += 0
        end
    end

    l = k_mais_l-k

    #Counting r: D2 is a matrix whose upper-right l-by-l block is diagonal,
    #   with the first r entries belonging to D2_til and the rest 1s.
    r = 0
    for i in 1:l
        if D2[i, k+i] != 1
            r += 1
        else
            r += 0
        end
    end

    bl = k
    br = l - r
    wl = m - k - r
    wr = p - l

    return bl,br,wl,wr,r
end

function permutation(wr,r,br, D2_julia)
    P = [zeros(r,wr)  I(r)        zeros(r,br);
         zeros(br,wr) zeros(br,r) I(br);
         I(wr)        zeros(wr,r) zeros(wr,br)]

    P_til = [zeros(wr,r) zeros(wr,br) I(wr);
             I(r)        zeros(r,br)  zeros(r,wr);
             zeros(br,r) I(br)        zeros(br,wr)]

    D2_our = P_til * D2_julia

    return P, D2_our
end


function op_switch(D1_til, D2_til, bl, br, wl, wr, r)
    # Cria novas matrizes com zeros
    D1_new = zeros(bl + r + wl, bl + r + wl)
    D2_new = zeros(wr + r + br, bl + r + wl)
  
    D1_new[1:bl, 1:bl] .= I(bl)
    D1_new[bl+1:bl+r, bl+1:bl+r] .= D1_til
    D1_new[bl+r+1:end, bl+r+1:end] .= I(wl)
  
    D2_new[wr+1:wr+r, bl+1:bl+r] .= D2_til
  
    return D1_new, D2_new
end


function aprox_by_angle(D1_til, D2_til, angle_lower_degrees::Float64, angle_upper_degrees::Float64, target_angle_degrees_for_keeping::Float64)
    if size(D1_til) != size(D2_til); error("D1_til e D2_til devem ter as mesmas dimensões"); end
    if !(0.0 <= angle_lower_degrees <= angle_upper_degrees <= 90.0); error("Limiares de ângulo inválidos. Deve ser 0 <= lower <= upper <= 90."); end
    if !(0.0 <= target_angle_degrees_for_keeping <= 90.0); error("Ângulo alvo inválido."); end
    
    r_val = size(D1_til, 1)
    if r_val == 0; return Diagonal(zeros(Float64,0)), Diagonal(zeros(Float64,0)); end

    D1_vec_aprox = zeros(Float64, r_val)
    D2_vec_aprox = zeros(Float64, r_val)

    target_angle_rad = deg2rad(target_angle_degrees_for_keeping)
    cos_target = cos(target_angle_rad)
    sin_target = sin(target_angle_rad)

    for i in 1:r_val
        c = D1_til[i,i]
        s = D2_til[i,i]
        angle_rad = atan(s, c) # Assumindo c,s >= 0 da GSVD, ângulo em [0, pi/2]
        angle_deg = rad2deg(angle_rad)
        
        #println("aprox_by_angle: Comp $i, c=$c, s=$s, ang_deg=$angle_deg")
        if angle_deg >= angle_lower_degrees && angle_deg <= angle_upper_degrees
            D1_vec_aprox[i] = cos_target
            D2_vec_aprox[i] = sin_target
        else
            D1_vec_aprox[i] = 1.0 # Para descarte, T_ii = 0
            D2_vec_aprox[i] = 0.0
        end
    end
    return Diagonal(D1_vec_aprox), Diagonal(D2_vec_aprox)
end

function op_aprox_by_angle(D1_til_gsvd::AbstractMatrix, D2_til_gsvd::AbstractMatrix, 
                           D1_structural_template::AbstractMatrix, D2_structural_template::AbstractMatrix, 
                           bl::Int, br::Int, wl::Int, wr::Int, r_val::Int, 
                           angle_lower_degrees::Float64, angle_upper_degrees::Float64, 
                           target_angle_degrees_for_keeping::Float64)

    # Criar cópias das matrizes estruturais para não modificar as originais
    D1_new = copy(D1_structural_template)
    D2_new = copy(D2_structural_template)

    if r_val > 0
        # Calcular os blocos _til aproximados
        D1_til_aprox, D2_til_aprox = aprox_by_angle(D1_til_gsvd, D2_til_gsvd, 
                                                    angle_lower_degrees, angle_upper_degrees, 
                                                    target_angle_degrees_for_keeping)

        # Verificar se as dimensões dos blocos _til aproximados são compatíveis
        # com as fatias onde serão inseridos.
        # D1_structural_template (e D1_new) tem o bloco _til em [bl+1:bl+r_val, bl+1:bl+r_val]
        # D2_structural_template (e D2_new) tem o bloco _til em [wr+1:wr+r_val, bl+1:bl+r_val]
        
        # Assegurar que D1_til_aprox e D2_til_aprox são r_val x r_val
        if size(D1_til_aprox) != (r_val, r_val) || size(D2_til_aprox) != (r_val, r_val)
            error("Dimensões de D1_til_aprox ou D2_til_aprox não são r_val x r_val como esperado.")
        end

        # Inserir os blocos _til aproximados nas posições corretas das cópias
        # Para D1_new, o bloco D1_til está em [bl+1 : bl+r_val, bl+1 : bl+r_val]
        if bl + r_val <= size(D1_new, 1) && bl + r_val <= size(D1_new, 2)
            D1_new[bl+1 : bl+r_val, bl+1 : bl+r_val] = D1_til_aprox
        elseif r_val > 0
             error("Dimensões de D1_new não comportam D1_til_aprox na posição esperada.")
        end

        # Para D2_new, o bloco D2_til está em [wr+1 : wr+r_val, bl+1 : bl+r_val]
        if wr + r_val <= size(D2_new, 1) && bl + r_val <= size(D2_new, 2)
            D2_new[wr+1 : wr+r_val, bl+1 : bl+r_val] = D2_til_aprox
        elseif r_val > 0
            error("Dimensões de D2_new não comportam D2_til_aprox na posição esperada.")
        end
    end
    # Se r_val == 0, D1_new e D2_new já são as cópias de D1_structural_template e D2_structural_template,
    # que devem estar corretas (sem blocos _til para modificar).

    return D1_new, D2_new
end



function run_test_aprox_by_angle_on_generated_data(A::AbstractMatrix, B::AbstractMatrix, test_case_name::String, angle_low::Float64, angle_high::Float64, target_angle::Float64)
 
    println("\n--- Executando Teste para: $test_case_name ---")
    println("Usando limiares de ângulo: [$angle_low, $angle_high] graus, Ângulo alvo para manter: $target_angle graus.")
  
    println("Matriz A (Entrada):"); display(A)
    println("Matriz B (Entrada):"); display(B)
  
    U, V_til, H_gsvd, D1_til_gsvd, D2_til_gsvd, D1_full_gsvd, D2_full_our_gsvd, bl_gsvd, br_gsvd, wl_gsvd, wr_gsvd, r_val_gsvd, P_perm_gsvd, V_orig_gsvd = Our_SVD(A, B)

    println("D1_full (Recuperado por Our_SVD):"); display(D1_full_gsvd)
    println("D2_full (Recuperado por Our_SVD):"); display(D2_full_our_gsvd)

    println("Parâmetros de bloco recuperados: bl=$bl_gsvd, br=$br_gsvd, wl=$wl_gsvd, wr=$wr_gsvd, r=$r_val_gsvd")

    D1_switch, D2_switch =  op_switch( D1_til_gsvd, D2_til_gsvd,bl_gsvd, br_gsvd, wl_gsvd, wr_gsvd, r_val_gsvd) # troca de cor e garante que D1 é inversível
    
    println("D1_switch (reconstrução -- troca de cor):"); display(D1_switch)
    println("D2_switch (reconstrução -- troca de cor):"); display(D2_switch)
  
    D1_new, D2_new = op_aprox_by_angle(D1_til_gsvd, D2_til_gsvd, D1_switch, D2_switch, 
                                       bl_gsvd, br_gsvd, wl_gsvd, wr_gsvd, r_val_gsvd, 
                                       angle_low, angle_high, target_angle)
    

    println("D1_new (após aproximação e reconstrução):"); display(D1_new)
    println("D2_new (após aproximação e reconstrução):"); display(D2_new)
    T_approx = D2_new * inv(D1_new)
    println("T_approx:"); display(T_approx)
  
    # Ajuste: usar somente as primeiras r colunas de U e V_til para compor Q_estimated,
    # de forma que as dimensões de T_approx (r x r) sejam compatíveis.
    r_cols = size(T_approx, 1)
    if size(V_til, 2) >= r_cols && size(U, 2) >= r_cols
        V_til_r = V_til[:, 1:r_cols]
        U_r = U[:, 1:r_cols]
        Q_estimated = V_til_r * T_approx * U_r'
        println("Q_estimated:"); display(Q_estimated)
  
        A_transformed = Q_estimated * A
        println("A_transformed:"); display(A_transformed)
  
        erro_frobenius = norm(B - A_transformed)
        @printf "Erro de Aproximação Global (Norma de Frobenius ||B - A_transformed||): %.4f\n" erro_frobenius
  
        println("\nComparação Coluna a Coluna (Norma da Diferença):")
        if size(A_transformed, 2) == size(B, 2)
            for j in 1:size(B, 2)
                norm_diff_col = norm(B[:,j] - A_transformed[:,j])
                @printf "Coluna %d: %.4f\n" j norm_diff_col
            end
        else
            println("Número de colunas de B e A_transformed não coincide.")
        end
    else
        println("Dimensões incompatíveis para truncamento em U/V_til. Pulando Q_estimated.")
    end
          # --- Comparação com Procrustes Absoluto ---
      Q_procrustes = absolute_procrustes(A, B)
      A_transformed_p = Q_procrustes * A
      erro_procrustes = norm(B - A_transformed_p)

      println("\n[Comparação de Métodos]")
      @printf "Erro Aproximação GSVD (Frobenius): %.6f\n" erro_frobenius
      @printf "Erro Aproximação Procrustes:      %.6f\n" erro_procrustes

      # Diferença entre transformações
      Q_diff = Q_estimated - Q_procrustes
      @printf "Norma da diferença entre Q_GSVD e Q_Procrustes: %.6f\n" norm(Q_diff)
  
    println("--- Fim do Teste: $test_case_name ---")
end

function generate_data_case(case_type::String, m_rows::Int, p_rows::Int, n_cols_H::Int, 
                            bl::Int, br::Int, wl::Int, wr::Int, r::Int, noise_level::Float64 = 0.01)
    println("\n--- Gerando Dados para: $case_type ---")
  
    dim_U_cols = bl + r + wl
    dim_V_cols = wr + r + br
    dim_H_rows = bl + r + br
  
    temp_U = randn(m_rows, dim_U_cols)
    U_ortho = Matrix(qr(temp_U).Q)
    if size(U_ortho,2) < dim_U_cols
        println("Advertência: m_rows < dim_U_cols, U pode não ter posto de coluna completo como desejado.")
    end
  
    temp_V = randn(p_rows, dim_V_cols)
    V_ortho = Matrix(qr(temp_V).Q)
    if size(V_ortho,2) < dim_V_cols
        println("Advertência: p_rows < dim_V_cols, V pode não ter posto de coluna completo como desejado.")
    end
  
    if dim_H_rows != n_cols_H
        println("Para H ser quadrada e invertível facilmente, defina n_cols_H = bl+r+br ($dim_H_rows).")
        if dim_H_rows <= n_cols_H
            H_matrix = randn(dim_H_rows, n_cols_H)
        else
            H_matrix = randn(dim_H_rows, n_cols_H) 
            println("Advertência: dim_H_rows > n_cols_H, H não terá posto de linha completo.")
        end
    else
         H_matrix = randn(dim_H_rows, n_cols_H)
         if rank(H_matrix) < dim_H_rows 
             H_matrix += Matrix(I,dim_H_rows,dim_H_rows)*0.1
         end
         if rank(H_matrix) < dim_H_rows 
             H_matrix = Matrix(I,dim_H_rows,dim_H_rows); 
             println("H forçada para Identidade.")
         end
    end
  
    D1_til_gen = zeros(r, r)
    D2_til_gen = zeros(r, r)
  
    if case_type == "Perfeito com Ruído"
        if r > 0
            D1_til_gen = Diagonal(fill(cosd(45), r))
            D2_til_gen = Diagonal(fill(sind(45), r))
        end
        A_gen, B_gen = Make_A_B(bl, br, wl, wr, r, D1_til_gen, D2_til_gen, H_matrix, U_ortho, V_ortho)
        A_final = A_gen + noise_level * randn(size(A_gen))
        B_final = B_gen + noise_level * randn(size(B_gen))
        println("D1_til usado na geração (Perfeito):"); display(D1_til_gen)
        println("D2_til usado na geração (Perfeito):"); display(D2_til_gen)
  
    elseif case_type == "Misto"
        if r >= 1
            D1_til_gen[1,1] = cosd(45); D2_til_gen[1,1] = sind(45)
        end
        if r >= 2
            D1_til_gen[2,2] = cosd(10); D2_til_gen[2,2] = sind(10)
        end
        for i in 3:r
            D1_til_gen[i,i] = cosd(45); D2_til_gen[i,i] = sind(45)
        end
  
        A_final, B_final = Make_A_B(bl, br, wl, wr, r, D1_til_gen, D2_til_gen, H_matrix, U_ortho, V_ortho)
        A_final += noise_level * randn(size(A_final))
        B_final += noise_level * randn(size(B_final))
        println("D1_til usado na geração (Misto):"); display(D1_til_gen)
        println("D2_til usado na geração (Misto):"); display(D2_til_gen)
    else
        error("Tipo de caso desconhecido: $case_type")
    end
    return A_final, B_final
end

function absolute_procrustes(A, B)
    U, _, Vt = svd(B * A')  # B Aᵀ e não Aᵀ B
    return U * Vt
end

function main_test_suite(A_input=nothing, B_input=nothing)
    angle_low_thresh = 30.0 
    angle_high_thresh = 60.0
    target_angle_keep = 45.0

    if A_input !== nothing && B_input !== nothing
        run_test_aprox_by_angle_on_generated_data(A_input, B_input, "Entrada Externa", 
                                                  angle_low_thresh, angle_high_thresh, target_angle_keep)
        return
    end

    # Código antigo continua aqui se nenhuma entrada for fornecida
    m_rows = 2  
    p_rows = 2 
    n_cols_H = 3   
    r_param = 1
    bl_param = 1
    br_param = 0
    n_cols_H_adjusted = bl_param + r_param + br_param
    if n_cols_H_adjusted == 0
        n_cols_H_adjusted = 1
    end
    if n_cols_H != n_cols_H_adjusted
        println("Ajustando n_cols_H para $n_cols_H_adjusted para H ser quadrada com D1/D2.")
        n_cols_H = n_cols_H_adjusted
    end
  
    wl_param = m_rows - bl_param - r_param 
    wr_param = p_rows - br_param - r_param
    if wl_param < 0; wl_param = 0; println("wl_param ajustado para 0"); end
    if wr_param < 0; wr_param = 0; println("wr_param ajustado para 0"); end

    println("Parâmetros de Geração: m=$m_rows, p=$p_rows, n_cols_AB=$n_cols_H")
    println("Blocos: bl=$bl_param, br=$br_param, wl=$wl_param, wr=$wr_param, r=$r_param")

    A_p, B_p = generate_data_case("Perfeito com Ruído", m_rows, p_rows, n_cols_H, 
                                  bl_param, br_param, wl_param, wr_param, r_param, 0.05)
    run_test_aprox_by_angle_on_generated_data(A_p, B_p, "Perfeito com Ruído", 
                                              angle_low_thresh, angle_high_thresh, target_angle_keep)

    A_m, B_m = generate_data_case("Misto", m_rows, p_rows, n_cols_H, 
                                  bl_param, br_param, wl_param, wr_param, r_param, 0.05)
    run_test_aprox_by_angle_on_generated_data(A_m, B_m, "Misto", 
                                              angle_low_thresh, angle_high_thresh, target_angle_keep)
end

#main_test_suite()

B = [0 1 0
     0 0 1]

A = [1 0 0
     0 1 0]

main_test_suite(A, B)
