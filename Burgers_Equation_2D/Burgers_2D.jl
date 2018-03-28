include("Thomas5DB.jl") # Algoritmo de Thomas, bloques tridiagonal

# Dominio mxm y nxn nodos desconocidos

mx = 50       # Es el numero de nodos en x
my = mx       # Es el numero de nodos en y
nx = mx-2      # Nodos desconocidos por lado en x
ny = my-2      # Nodos desconocidos por lado en y
x0 = 0
y0 = 0
xf = 1
yf = 1
dx = (xf-x0)/(mx-1)
dy = (yf-y0)/(my-1)

R = 22
miu = 1/R
dt = 0.0001  
iter=1000

# Se definen el tipo de condiciones de frontera
# Neuman du/dx=0; dv/dy=0
# Si p=1 entonces es condicion Neuman
# Si p=0 entonces es condicion Dirichlet

p1 = 0   # Solo aplica para u
p3 = 0   # Solo aplica para u

p2 = 0   # Solo aplica para v
p4 = 0   # Solo aplica para v


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

Mu0=zeros(Float64,mx,my)
Mv0=zeros(Float64,mx,my)
Mh0=zeros(Float64,mx,my)

for i in 1:mx
   for j in 1:my
     x=VecX[i]
     y=VecY[j]
     exponente = R*(-4*x+4*y)/32
     Mu0[i,j] = (3/4)-1/(4*(1+e^(exponente)))
     Mv0[i,j] = (3/4)+1/(4*(1+e^(exponente)))
   end
end

println("Condiciones iniciales creadas")

# -------- VecU y VecV reorganizan en un vector para crear las constantes
VecU = zeros(Float64,nx*ny)
VecV = zeros(Float64,nx*ny)

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

for i in 1:mx
for j in 1:my
   VecXPlot[j+(i-1)*my] = VecX[i]
   VecYPlot[j+(i-1)*my] = VecY[j]
   MsolU[j+(i-1)*my, 1] = Mu0[i,j]
   MsolV[j+(i-1)*my, 1] = Mv0[i,j]
end
end


for i in 1:nx
for j in 1:ny
   VecU[j+(i-1)*ny]=Mu0[i+1,j+1] # mxm
   VecV[j+(i-1)*ny]=Mv0[i+1,j+1] # mxm
end
end

# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -------------- Este ciclo se repite cada paso de tiempo ---
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------

for k in 2:iter+1 # Inicio del gran ciclo

t=Vect[k]
println(t)

# Diagonales de la matriz que se resolvera en cada iteracion
Mx1 = zeros(Float64,nx*ny-ny) # Diagonal inferior mas alejada
Mx2 = zeros(Float64,nx*ny-1)  # Diagonal inferior segunda
Mx3 = zeros(Float64,nx*ny)    # Diagonal principal
Mx4 = zeros(Float64,nx*ny-1)  # Diagonal superior segunda
Mx5 = zeros(Float64,nx*ny-ny) # Diagonal superior mas alejada

My1 = zeros(Float64,nx*ny-ny) # Diagonal inferior mas alejada
My2 = zeros(Float64,nx*ny-1)  # Diagonal inferior segunda
My3 = zeros(Float64,nx*ny)    # Diagonal principal
My4 = zeros(Float64,nx*ny-1)  # Diagonal superior segunda
My5 = zeros(Float64,nx*ny-ny) # Diagonal superior mas alejada

   for i in 1:mx
     x=VecX[i]

     exp2 = R*(-4*x + 4*yf -t)/32 # Pared 2
     exp4 = R*(-4*x + 4*y0 -t)/32 # Pared 4

     CC2[i,1] = (3/4)-1/(4*(1+e^(exp2)))
     CC4[i,1] = (3/4)-1/(4*(1+e^(exp4)))

     CC2[i,2] = (3/4)+1/(4*(1+e^(exp2)))
     CC4[i,2] = (3/4)+1/(4*(1+e^(exp4)))
   end

   for i in 1:my
     y=VecY[i]

     exp1 = R*(-4*x0+ 4*y  -t)/32 # Pared 1
     exp3 = R*(-4*xf+ 4*y  -t)/32 # Pared 3

     CC1[i,1] = (3/4)-1/(4*(1+e^(exp1)))
     CC3[i,1] = (3/4)-1/(4*(1+e^(exp3)))

     CC1[i,2] = (3/4)+1/(4*(1+e^(exp1)))
     CC3[i,2] = (3/4)+1/(4*(1+e^(exp3)))
   end

for i in 1:nx*ny # Se escriben las diag. prin.
     B = -2*miu/(dx^2) -2*miu/(dy^2) -1/dt # Cte diagonal pr.

     Mx3[i] = B
     My3[i] = B
end


for i in 1:nx*ny-1  # Se escriben la diag. segun. sup.
     u_ij = VecU[i] # Termino no lineal U
     v_ij = VecV[i] # Termino no lineal V

     E = miu/(dy^2) - v_ij/(2*dy)

     Mx4[i]=E
     My4[i]=E
end

for i in 1:nx*ny-1 # Se escriben la diag. segun. inf.
     u_ij = VecU[i+1] # Termino no lineal U
     v_ij = VecV[i+1] # Termino no lineal V

     D = miu/(dy^2) + v_ij/(2*dy)

     Mx2[i] = D
     My2[i] = D
end

#for i in 1:nx-1 # Ceros en las diag. secundarias
 #    Mx4[ny*i,i*ny+1] = 0.00000000  # Superior
 #    Mx2[i*ny+1,ny*i] = 0.00000000  # Inferior
 #    My4[ny*i,i*ny+1] = 0.00000000  # Superior
 #    My2[i*ny+1,ny*i] = 0.00000000  # Inferior
#end

for i in 1:nx*ny-ny  # Se escriben la diag. mas sup.
     u_ij = VecU[i] # Termino no lineal U
     v_ij = VecV[i] # Termino no lineal V

     C = miu/(dx^2) - u_ij/(2*dx)

     Mx5[i]=C
     My5[i]=C
end

for i in 1:nx*ny-ny  # Se escriben la diag. mas inf.
     u_ij = VecU[i+ny] # Termino no lineal U
     v_ij = VecV[i+ny] # Termino no lineal V

     A = miu/(dx^2) + u_ij/(2*dx)

     Mx1[i] = A
     My1[i] = A
end

# -------- Vectores conocidos se reinician
VecBx = zeros(Float64,nx*ny)
VecBy = zeros(Float64,nx*ny)

VecBx = (-1/dt)*VecU  # Vector B conocido en X
VecBy = (-1/dt)*VecV  # Vector B conocido en Y

# ---------------------------------------------------------------
# -------------- Se escriben las condiciones Dirichet -----------
# ---------------------------------------------------------------

for i in 1:ny # Se escriben las CC de la pared 1
     u_ij = CC1[i+1,1]   # Nodos pared 1 en x
     v_ij = CC1[i+1,2]   # Nodos pared 1 en y

     u_nl= VecU[i] # Termino no lineal U
     v_nl= VecV[i] # Termino no lineal V

     A = miu/(dx^2) + u_nl/(2*dx)

     VecBx[i] = VecBx[i] - A*u_ij*(1-p1)
     VecBy[i] = VecBy[i] - A*v_ij
end

for i in 1:ny # Se escriben las CC de la pared 3
     u_ij = CC3[i+1,1]   # Nodos pared 3 en x
     v_ij = CC3[i+1,2]   # Nodos pared 3 en y

     u_nl= VecU[i+nx*ny-ny] # Termino no lineal U
     v_nl= VecV[i+nx*ny-ny] # Termino no lineal V

     C = miu/(dx^2) - u_nl/(2*dx)

     VecBx[i+nx*ny-ny] = VecBx[i+nx*ny-ny] - C*u_ij*(1-p3)
     VecBy[i+nx*ny-ny] = VecBy[i+nx*ny-ny] - C*v_ij
end

for i in 1:nx # Se escriben las CC de la pared 2
     u_ij = CC2[i+1,1]   # Nodos pared 2 en x
     v_ij = CC2[i+1,2]   # Nodos pared 2 en y

     u_nl= VecU[i*ny] # Termino no lineal U
     v_nl= VecV[i*ny] # Termino no lineal V

     E = miu/(dy^2) - v_nl/(2*dy)

     VecBx[i*ny] = VecBx[i*ny] - E*u_ij
     VecBy[i*ny] = VecBy[i*ny] - E*v_ij*(1-p2)
end

for i in 0:nx-1 # Se escriben las CC de la pared 4
     u_ij = CC4[i+2,1]   # Nodos pared 4 en x
     v_ij = CC4[i+2,2]   # Nodos pared 4 en y

     u_nl= VecU[i*ny+1] # Termino no lineal U
     v_nl= VecV[i*ny+1] # Termino no lineal V

     D = miu/(dy^2) + v_nl/(2*dy)

     VecBx[i*ny+1] = VecBx[i*ny+1] - D*u_ij
     VecBy[i*ny+1] = VecBy[i*ny+1] - D*v_ij*(1-p4)
end

#---------------------------------------------------------------
#-----------------Solucion de los dos sistemas -----------------
#---------------------------------------------------------------

Thomas5G(Mx1,Mx2,Mx3,Mx4,Mx5,VecBx,VecU,nx)
Thomas5G(My1,My2,My3,My4,My5,VecBy,VecV,nx)

    for i in 1:nx
    for j in 1:ny # Se guardan los datos en Msol
      MsolU[j+1+my*(i),k] = VecU[(i-1)*ny+j]
      MsolV[j+1+my*(i),k] = VecV[(i-1)*ny+j]
    end
    end

# Se agregan las condiciones de contorno

    for i in 1:mx
      MsolU[i*my,k] = CC2[i,1] 
      MsolU[(i-1)*my+1,k] = CC4[i,1] 

      MsolV[i*my,k] = CC2[i,2] 
      MsolV[(i-1)*my+1,k] = CC4[i,2] 
    end

    for i in 1:my
      MsolU[i,k]   = CC1[i,1] 
      MsolU[my*(mx-1)+i,k] = CC3[i,1] 

      MsolV[i,k]   = CC1[i,2] 
      MsolV[my*(mx-1)+i,k] = CC3[i,2] 
    end


# Datos que se tomaran en la siguiente iteracion
#for i in 1:mx
#for j in 1:my
  # Mu0[i,j]=MsolU[j+(i-1)*my , k]
 #  Mv0[i,j]=MsolV[j+(i-1)*my , k]
#end
#end

end # Final del ciclo Principal

println("Errores")

# --------------------------------------------------------------
# ------------------ Se calcula la solucion exacta -------------
# --------------------------------------------------------------

    M_exa_x = zeros(Float64,mx*my) # Matriz de soluciones en X(t)
    M_exa_y = zeros(Float64,mx*my) # Matriz de soluciones en Y(t)
    Mx0=zeros(Float64,mx,my)
    My0=zeros(Float64,mx,my)
    vt=zeros(Float64,5)
    Err_x = zeros(Float64,mx*my) # Vector de errores relativos en x
    Err_y = zeros(Float64,mx*my) # Vector de errores relativos en y

vt=[10 100 1000 5000 10000 20000]

for k in 1:1

     t = vt[k]*dt # Tiempo para la solucion analitica
     println(t)

    for i in 1:mx
       for j in 1:my
        x=VecX[i]
        y=VecY[j]
        exponente = R*(-4*x+4*y-t)/32
        Mx0[i,j] = (3/4)-1/(4*(1+e^(exponente)))
        My0[i,j] = (3/4)+1/(4*(1+e^(exponente)))
       end
    end


    for i in 1:mx
    for j in 1:my
       M_exa_x[j+(i-1)*my] = Mx0[i,j]
       M_exa_y[j+(i-1)*my] = My0[i,j]
    end
    end

    

    for i in 1:mx*my
     Err_x[i] = (MsolU[i,vt[k]+1] - M_exa_x[i])
    end
println(maximum(Err_x))
end


println(dx)
# --------------------------------------------------------------
# --------------------- Se guardan los datos -------------------
# --------------------------------------------------------------


open("VelocidadU.dat","w") do io
     writedlm(io,[VecXPlot VecYPlot MsolU])
end

open("VelocidadV.dat","w") do io
     writedlm(io,[VecXPlot VecYPlot MsolV])
end


println("Terminado")


