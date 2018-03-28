El archivo "Burgers_2D.jl" toma un dominio rectangular, condiciones Dirichlet o Neumann para la velocidad en las paredes. 
Se generan dos archivos con las velocidades U y V en cada instante de tiempo.
Se tiene un sistema tridiagonal por bloques (Block-tridiagonal matrix), el cual se resuelve con la 
función definida en el archivo "Thomas5DB.jl"
El archivo "B2D.gp" da las condiciones al archivo "B2D.plt" para crear imágenes de cada instante de tiempo 
y guardarlas en una carpeta "animation".
