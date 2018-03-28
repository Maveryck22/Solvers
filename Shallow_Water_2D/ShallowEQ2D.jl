

# Dominio mxm y nxn nodos desconocidos

mx = 200       # Es el numero de nodos en x
my = 200       # Es el numero de nodos en y
nx = mx-2      # Nodos desconocidos por lado en x
ny = my-2      # Nodos desconocidos por lado en y
x0 = 0
y0 = 0
xf = 0.6
yf = 0.6
dx = (xf-x0)/(mx-1)
dy = (yf-y0)/(my-1)

Kx = 0.0002
Ky = 0.0002
g = 10
dt = 0.001  
iter=0

# Vectores de coordenadas

VecX = zeros(Float64,mx)
VecY = zeros(Float64,my)
Vect = zeros(Float64,iter+1)

for i in 1:mx
  VecX[i]= x0 + dx*(i-1)
end

for i in 1:my
  VecY[i]= y0 + dy*(i-1)
end

for i in 1:iter+1
  Vect[i]=dt*(i-1)
end


# Matrices con las condiciones iniciales

Mu0= zeros(Float64,mx,my)
Mv0= zeros(Float64,mx,my)
Mz0= zeros(Float64,mx,my)
Mh = zeros(Float64,mx,my)

for i in 1:mx
   for j in 1:my
     x = VecX[i]
     y = VecY[j]
      Soliton1 = (3.2)*e^(-80*(x-0.2)^2-80*(y-0.2)^2)
      Soliton2 = (1.8)*e^(-100*(x-0.3)^2-100*(y-0.3)^2)
     Mz0[i,j] = 1 + Soliton2 #+ Soliton2
      Mh[i,j] = 0
   end
end

println("Condiciones iniciales creadas")

# -------- VecU y VecV reorganizan en un vector para crear las constantes
VecU = zeros(Float64,nx*ny)
VecV = zeros(Float64,nx*ny)
VecZ = zeros(Float64,nx*ny)
VecH2= zeros(Float64,nx*ny)
VecH = zeros(Float64,mx*my)

CC1 = zeros(Float64,my,2) # Cond. cont. pared 1
CC2 = zeros(Float64,mx,2) # Cond. cont. pared 2
CC3 = zeros(Float64,my,2) # Cond. cont. pared 3
CC4 = zeros(Float64,mx,2) # Cond. cont. pared 4

# ---------  Vectores para graficar

VecXPlot = zeros(Float64,mx*my)
VecYPlot = zeros(Float64,mx*my)

MsolU=zeros(Float64,mx*my,iter+1) # Matriz con las soluciones de U(t)
MsolV=zeros(Float64,mx*my,iter+1) # Matriz con las soluciones de V(t)
MsolZ=zeros(Float64,mx*my,iter+1) # Matriz con las soluciones de z(t)

VecUSol = zeros(Float64,nx*ny)
VecVSol = zeros(Float64,nx*ny)
VecZSol = zeros(Float64,nx*ny)

for i in 1:mx
for j in 1:my
   VecXPlot[j+(i-1)*my] = VecX[i]
   VecYPlot[j+(i-1)*my] = VecY[j]
   MsolZ[j+(i-1)*my, 1] = Mz0[i,j]
end
end

for i in 1:nx
for j in 1:ny
   VecU[j+(i-1)*ny] = Mu0[i+1,j+1] # mxm
   VecV[j+(i-1)*ny] = Mv0[i+1,j+1] # mxm
   VecZ[j+(i-1)*ny] = Mz0[i+1,j+1] # mxm
  VecH2[j+(i-1)*ny] = Mh[i+1,j+1] # mxm
end
end

for i in 1:mx
for j in 1:my
   VecH[j+(i-1)*my] =  Mh[i,j] # mxm
end
end


MatA = zeros(Float64,nx*ny,nx*ny)
MatB = zeros(Float64,nx*ny,nx*ny)
MatC = zeros(Float64,nx*ny,nx*ny)
MatD = zeros(Float64,nx*ny,nx*ny)
MatE = zeros(Float64,nx*ny,nx*ny)
MatH = zeros(Float64,nx*ny,nx*ny)
MatI = zeros(Float64,nx*ny,nx*ny)

VecJ = zeros(Float64,nx*ny)
VecF = zeros(Float64,nx*ny)
VecG = zeros(Float64,nx*ny)

# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -------------- Este ciclo se repite cada paso de tiempo ---
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------

for k in 2:iter+1 # Inicio del ciclo PRINCIPAL

MatA = zeros(Float64,nx*ny,nx*ny)
MatB = zeros(Float64,nx*ny,nx*ny)
MatC = zeros(Float64,nx*ny,nx*ny)
MatD = zeros(Float64,nx*ny,nx*ny)
MatE = zeros(Float64,nx*ny,nx*ny)
MatH = zeros(Float64,nx*ny,nx*ny)
MatI = zeros(Float64,nx*ny,nx*ny)

VecJ = zeros(Float64,nx*ny)
VecF = zeros(Float64,nx*ny)
VecG = zeros(Float64,nx*ny)

t = Vect[k]
println(t)

   for i in 1:mx
     x=VecX[i]

     CC2[i,1] = 0
     CC4[i,1] = 0

     CC2[i,2] = 0
     CC4[i,2] = 0
   end

   for i in 1:my
     y=VecY[i]

     CC1[i,1] = 0
     CC3[i,1] = 0

     CC1[i,2] = 0
     CC3[i,2] = 0
   end

for i in 1:nx*ny # Se escriben las diag. prin.
     h11 = VecH[i]
     h21 = VecH[i+my]
     h12 = VecH[i+1]

     MatE[i,i] = 0 #(h21-h11)/(2*dx)   # E
     MatH[i,i] = 0 #(h12-h11)/(2*dy)   # H
     MatI[i,i] = 1/dt + 2*Kx/(dx^2) + 2*Ky/(dy^2)   # B
     MatA[i,i] = 1/dt      # G1
     MatC[i,i] = 1/dt      # G2
end

for i in 1:nx*ny-1  # Se escriben la diag. segun. sup.
     u_ij = VecU[i] # Termino no lineal U
     v_ij = VecV[i] # Termino no lineal V
     h_ij = VecH2[i]
     z_ij = VecZ[i]

     MatH[i,i+1] = (h_ij+z_ij)/(2*dy)     # I
     MatI[i,i+1] = v_ij/(2*dy) -Ky/(dy^2) # K
     MatA[i,i+1] = v_ij/(2*dy)            # F1
     MatC[i,i+1] = v_ij/(2*dy)            # D2
     MatD[i,i+1] = g/(2*dy)               # B2
end

for i in 1:nx*ny-1 # Se escriben la diag. segun. inf.
     u_ij = VecU[i+1] # Termino no lineal U
     v_ij = VecV[i+1] # Termino no lineal V
     h_ij = VecH2[i+1]
     z_ij = VecZ[i+1]

     MatH[i+1,i] = -(h_ij+z_ij)/(2*dy)     # G
     MatI[i+1,i] = -v_ij/(2*dy) -Ky/(dy^2) # J
     MatA[i+1,i] = -v_ij/(2*dy)            # E1
     MatC[i+1,i] = -v_ij/(2*dy)            # C2
     MatD[i+1,i] = -g/(2*dy)               # A2
end

for i in 1:nx-1 # Ceros en las diag. secundarias
     MatH[ny*i,i*ny+1] = 0.00000000  # Superior
     MatH[i*ny+1,ny*i] = 0.00000000  # Inferior
     MatI[ny*i,i*ny+1] = 0.00000000  # Superior
     MatI[i*ny+1,ny*i] = 0.00000000  # Inferior
     MatA[ny*i,i*ny+1] = 0.00000000  # Superior
     MatA[i*ny+1,ny*i] = 0.00000000  # Inferior
     MatC[ny*i,i*ny+1] = 0.00000000  # Superior
     MatC[i*ny+1,ny*i] = 0.00000000  # Inferior
     MatD[ny*i,i*ny+1] = 0.00000000  # Superior
     MatD[i*ny+1,ny*i] = 0.00000000  # Inferior
end

for i in 1:nx*ny-ny  # Se escriben la diag. mas sup.
     u_ij = VecU[i] # Termino no lineal U
     v_ij = VecV[i] # Termino no lineal V
     h_ij = VecH2[i]
     z_ij = VecZ[i]

     MatE[i,i+ny] = (h_ij+z_ij)/(2*dx)      # F
     MatI[i,i+ny] = u_ij/(2*dx) - Kx/(dx^2) # C
     MatA[i,i+ny] = u_ij/(2*dx)             # D1
     MatB[i,i+ny] = g/(2*dx)                # B1
     MatC[i,i+ny] = u_ij/(2*dx)             # F2
end

for i in 1:nx*ny-ny  # Se escriben la diag. mas inf.
     u_ij = VecU[i+ny] # Termino no lineal U
     v_ij = VecV[i+ny] # Termino no lineal V
     h_ij = VecH2[i+ny]
     z_ij = VecZ[i+ny]

     MatE[i+ny,i] = -(h_ij+z_ij)/(2*dx)      # D
     MatI[i+ny,i] = -u_ij/(2*dx) - Kx/(dx^2) # A
     MatA[i+ny,i] = -u_ij/(2*dx)             # C1
     MatB[i+ny,i] = -g/(2*dx)                # A1
     MatC[i+ny,i] = -u_ij/(2*dx)             # E2
end

# -------- Vectores conocidos se reinician

VecJ = (1/dt)*VecZ  # Vector B conocido en X
VecF = (1/dt)*VecU  # Vector B conocido en Y
VecG = (1/dt)*VecV  # Vector B conocido en Y

# ---------------------------------------------------------------
# -------------- Se escriben las condiciones Dirichet -----------
# ---------------------------------------------------------------

for i in 1:ny # Se escriben las CC de la pared 1
     u_ij = CC1[i+1,1]   # Nodos pared 1 en x
     v_ij = CC1[i+1,2]   # Nodos pared 1 en y

     u_nl= VecU[i] # Termino no lineal U
     v_nl= VecV[i] # Termino no lineal V

     VecF[i] = VecF[i] 
     VecG[i] = VecG[i] 
end

for i in 1:ny # Se escriben las CC de la pared 3
     u_ij = CC3[i+1,1]   # Nodos pared 3 en x
     v_ij = CC3[i+1,2]   # Nodos pared 3 en y

     u_nl= VecU[i+nx*ny-ny] # Termino no lineal U
     v_nl= VecV[i+nx*ny-ny] # Termino no lineal V

     VecF[i+nx*ny-ny] = VecF[i+nx*ny-ny] 
     VecG[i+nx*ny-ny] = VecG[i+nx*ny-ny] 
end

for i in 1:nx # Se escriben las CC de la pared 2
     u_ij = CC2[i+1,1]   # Nodos pared 2 en x
     v_ij = CC2[i+1,2]   # Nodos pared 2 en y

     u_nl= VecU[i*ny] # Termino no lineal U
     v_nl= VecV[i*ny] # Termino no lineal V

     VecF[i*ny] = VecF[i*ny] 
     VecG[i*ny] = VecG[i*ny] 
end

for i in 0:nx-1 # Se escriben las CC de la pared 4
     u_ij = CC4[i+2,1]   # Nodos pared 4 en x
     v_ij = CC4[i+2,2]   # Nodos pared 4 en y

     u_nl= VecU[i*ny+1] # Termino no lineal U
     v_nl= VecV[i*ny+1] # Termino no lineal V

     VecF[i*ny+1] = VecF[i*ny+1] 
     VecG[i*ny+1] = VecG[i*ny+1] 
end

# ---------------------------------------------------------------
# ---------------Se escriben las condiciones Neuman -------------
# ---------------------------------------------------------------

for i in 1:ny  # Para las paredes 1 y 3
     i2 = nx*ny-ny+i

     u1 = VecU[i]  # Termino no lineal U cercanos a la pared 1
     u3 = VecU[i2] # Termino no lineal U cercanos a la pared 3

     A = -Kx/(dx^2) - u1/(2*dx)
     C = -Kx/(dx^2) + u3/(2*dx)

     MatI[i,i]    = MatI[i,i]  + 4*A/3  # En la diagonal principal
     MatI[i,i+ny] = MatI[i,i+ny] - A/3  # En la diagonal mas superior

     MatI[i2,i2]    = MatI[i2,i2]  + 4*C/3 # En la diag. principal
     MatI[i2,i2-ny] = MatI[i2,i2-ny] - C/3 # En la diag. mas inferior

     A1 = -g/(2*dx)
     B1 =  g/(2*dx)

     MatB[i,i]    = MatB[i,i]  + 4*A1/3 # En la diag. principal
     MatB[i,i+ny] = MatB[i,i+ny] - A1/3  # En la diagonal mas superior

     MatB[i2,i2]    = MatB[i2,i2]  + 4*B1/3 # En la diag. principal
     MatB[i2,i2-ny] = MatB[i2,i2-ny] - B1/3 # En la diag. mas inferior
end

for i in 1:nx  # Para las paredes 2 y 4
     i4 = 1+(i-1)*ny
     i2 = i*ny

     v4 = VecV[i4]  # Termino no lineal v pared 4
     v2 = VecV[i2]  # Termino no lineal v pared 2

     J  = -Ky/(dy^2) - v4/(2*dy)
     K1 = -Ky/(dy^2) + v2/(2*dy)

     MatI[i4,i4]   = MatI[i4,i4] + 4*J/3  # En la diagonal principal
     MatI[i4,i4+1] = MatI[i4,i4+1] - J/3  # En la diagonal 2da sup.

     MatI[i2,i2]   = MatI[i2,i2] + 4*K1/3  # En la diagonal principal
     MatI[i2,i2-1] = MatI[i2,i2-1] - K1/3  # En la diagonal 2da inf.

     A2 = -g/(2*dy)
     B2 =  g/(2*dy)

     MatD[i4,i4]   = MatD[i4,i4] + 4*A2/3  # En la diagonal principal
     MatD[i4,i4+1] = MatD[i4,i4+1] - A2/3  # En la diagonal 2da sup.

     MatD[i2,i2]   = MatD[i2,i2] + 4*B2/3  # En la diagonal principal
     MatD[i2,i2-1] = MatD[i2,i2-1] - B2/3  # En la diagonal 2da inf.
end

#---------------------------------------------------------------
#-----------------Solucion de los dos sistemas -----------------
#---------------------------------------------------------------

A_1 = inv(MatA)
C_1 = inv(MatC)

M1 = VecJ - MatE*(A_1*VecF) - MatH*(C_1*VecG) 
M2 = MatI - MatE*(A_1*MatB) - MatH*(C_1*MatD) 

VecZ = inv(M2)*M1
VecU = A_1*(VecF - MatB*VecZ)
VecV = C_1*(VecG - MatD*VecZ)

    for i in 1:nx
    for j in 1:ny # Se guardan los datos en Msol
      MsolU[j+1+my*(i),k] = VecU[(i-1)*ny+j]
      MsolV[j+1+my*(i),k] = VecV[(i-1)*ny+j]
      MsolZ[j+1+my*(i),k] = VecZ[(i-1)*ny+j]
    end
    end

# Se agregan las CC para la velocidad (Dirichlet)

    for i in 1:mx
      MsolU[i*my,k] = 0#CC2[i,1] 
      MsolU[(i-1)*my+1,k] = 0#CC4[i,1] 

      MsolV[i*my,k] = 0#CC2[i,2] 
      MsolV[(i-1)*my+1,k] = 0#CC4[i,2] 
    end

    for i in 1:my
      MsolU[i,k]   = 0#CC1[i,1] 
      MsolU[my*(mx-1)+i,k] = 0#CC3[i,1] 

      MsolV[i,k]   = 0#CC1[i,2] 
      MsolV[my*(mx-1)+i,k] = 0#CC3[i,2] 
    end

# Se agregan las CC para la altura (Neumann)

    for i in 2:ny+1 # Para las paredes 1 y 3
      l = mx*my-my+i

      MsolZ[i,k] = (4/3)*MsolZ[i+my,k]-(1/3)*MsolZ[i+2*my,k]
      MsolZ[l,k] = (4/3)*MsolZ[l-my,k]-(1/3)*MsolZ[l-2*my,k]
    end

    for i in 1:nx+2 # Para las paredes 2 y 4
       l = (i-1)*my+1
  
      MsolZ[i*my,k] = (4/3)*MsolZ[i*my-1,k]-(1/3)*MsolZ[i*my-2,k]
      MsolZ[l,k] =  (4/3)*MsolZ[l+1,k]-(1/3)*MsolZ[l+2,k]
    end

end # Final del ciclo Principal


# --------------------------------------------------------------
# --------------------- Se guardan los datos -------------------
# --------------------------------------------------------------


open("Alturas.dat","w") do io
     writedlm(io,[VecXPlot VecYPlot MsolZ])
end

open("VelocidadU.dat","w") do io
     writedlm(io,[VecXPlot VecYPlot MsolU])
end

open("VelocidadV.dat","w") do io
     writedlm(io,[VecXPlot VecYPlot MsolV])
end


println("Terminado")



