# RBD

Modelo para calcular trayectoria de cuerpo rígido 6 DOF utilizando cuaterniones para la representación matemática

ventaja del método -> con UNA sola corrida podemos obtener todos los coeficientes

Para correr por primera vez, hacer un pull y crear una carpeta Resultados (en ~/RBD/Resultados) donde el prog escribirá archivos de coeficientes, fuerzas. Además los procesará para guardar en Forces_coef_proc.dat un coeficiente por cada paso de tiempo, ya que el integrador temporal tiene un paso adaptativo y para un mismo paso de tiempo calcula los coeficientes, entonces con la rutina proc.py se reacomoda los resultados.