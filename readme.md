# Implementación de Filtro de Kalman Extendido (EKF) en C

Esta librería está en formato de componente para su uso con ESP-IDF.

Esta librería proporciona una implementación de un Filtro de Kalman Extendido (EKF) en C, optimizada para la fusión de sensores en sistemas de navegación y control. El EKF es un algoritmo recursivo que estima el estado de un sistema dinámico a partir de una serie de mediciones incompletas o ruidosas. 

Esta rama, incorpora corrección independiente de heading usando medidas de magnetómetro y un vector de estados extendido inncluyendo el bias del gisoscopio.
