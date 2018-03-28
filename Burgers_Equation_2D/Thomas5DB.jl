function Thomas5G(D1,D2,D3,D4,D5,VB,VX,nx)

# D1: Diagonal 1 inferior
# D2: Diagonal 2 inferior
# D3: Diagonal principal
# D4: Diagonal 1 superior
# D5: Diagonal 2 superior
# VB: Vector derecho
# VX: Vector solucion
# Los vectores originales no se modifican

  n2 = length(D1)  # Diagonales mas lejanas
  n = length(D3)  # Diagonal principal

  ny = n-n2

# Matrices donde se guardan los bloques-matrices
MC = zeros(Float64,(nx-1)*ny,ny) # Bloques superiores
VD = zeros(Float64,n)       # vectores derechos por bloques

# Matrices auxiliares
A = zeros(Float64,ny,ny)
B = zeros(Float64,ny,ny)
C = zeros(Float64,ny,ny)
D = zeros(Float64,ny)
X = zeros(Float64,ny)
IB= zeros(Float64,ny,ny)
IC= zeros(Float64,ny,ny)
ID= zeros(Float64,ny)

# Primera ecuación, bloque 1
   for i in 1:ny
      A[i,i] = D1[i]  # Bloques inferiores
      B[i,i] = D3[i]  # Bloques diagonales
      C[i,i] = D5[i]  # Bloques superiores
      D[i]   = VB[i]  # Bloques de vectores B
    end
    for i in 1:ny-1
       B[i+1,i] = D2[i]
       B[i,i+1] = D4[i]
    end
   IC = inv(B)*C
   ID = inv(B)*D

# Se guardan C1' y D1'
    for i in 1:ny 
    for j in 1:ny 
     MC[i,j] = IC[i,j]
    end
     VD[i] = ID[i] 
    end

# Se inicia el barrido hasta nx-1 bloques

  for k in 1:nx-1 # Forward sweep

        for i in 1:ny
           D[i]   = VB[i+k*ny]
           B[i,i] = D3[i+k*ny]
           A[i,i] = D1[i+(k-1)*ny]
         end

           if k<=nx-2
            for i in 1:ny 
             C[i,i] = D5[i+k*ny]
            end
           end

         for i in 1:ny-1
            B[i+1,i] = D2[i+k*ny]
            B[i,i+1] = D4[i+k*ny]
         end

         ID = inv(B-A*IC)*(D-A*ID)
         IC = inv(B-A*IC)*C

       # Se guardan Ci' y Di'
          for i in 1:ny 
            VD[i+k*ny] = ID[i] 
          end

          for i in 1:ny 
          for j in 1:ny 
            m2 = i+k*ny
            if m2<= n-ny
            MC[m2,j] = IC[i,j]
            end
          end
          end
   end

  # Se reemplaza el ultimo elemento
   X=ID
   for i in ny:-1:1
      VX[i+n-ny] = ID[i]
   end

     for k in (nx-1):-1:1 # Back substitution
       # Se leen Ci' y Di'
          for i in 1:ny
          for j in 1:ny 
            C[i,j] = MC[i+(k-1)*ny,j]
          end
            D[i] = VD[i+(k-1)*ny]
          end

        X = D-C*X

        for i in ny:-1:1 
         VX[i+(k-1)*ny] = X[i]
        end
     end

     return VX
end # Fin de la funcion

