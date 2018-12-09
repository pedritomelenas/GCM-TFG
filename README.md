# GCM-TFG Galaxy components modeling

Este repositorio contiene el código implementado para el TFM de Marina Estévez Almenzar, dirigido por
[Óscar Sánchez Romero](https://www.ugr.es/~ossanche) [@oscarsanchezromero](https://github.com/oscarsanchezromero) y [Pedro A. García Sánchez](https://wwww.ugr.es/~pedro) [@pedritomelenas](https://github.com/pedritomelenas).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pedritomelenas/GCM-TFG/master)

## Resumen

En este trabajo se presenta y desarrolla un problema de optimización presente actualmente en el campo de la Astrofísica. 
Dada una galaxia, se plantea el ajuste de su curva de rotación con el fin de determinar la cantidad de materia oscura que 
hay en dicha galaxia, con la mejor exactitud posible. Desde el punto de vista matemático, se nos presenta un problema de 
minimización de un funcional que, inicialmente, depende de más de un parámetro, y se expresa como un problema de mínimos 
cuadrados ponderado. Mediante un proceso de optimización anidada, y en virtud de las buenas propiedades matemáticas de uno 
de estos parámetros, conseguiremos reducir el espacio paramétrico, facilitando así el problema. Desde el punto de vista 
informático, será necesario desarrollar algoritmos que se adapten al problema. Se nos plantea la búsqueda del mínimo de 
una función, sin tener asegurada su existencia, en un espacio que inicialmente es infinito. En base a las propiedades 
del problema habrá que acotar convenientemente el espacio y buscar el mínimo en el mismo de forma eficiente.
