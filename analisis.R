
# libraries ---------------------------------------------------------------

library(tidyverse)
library(mixexp)
library(openxlsx)
library(desirability)
library(car)
library(gtools)
library(gridExtra)


# Funciones customizadas --------------------------------------------------

source("custom_MixturePlot.R")

# read data ---------------------------------------------------------------

datos <- openxlsx::read.xlsx("datos/datos_exfoliante.xlsx", cols = 2:9)

variables <- names(datos[, 1:6])


# analyze data  -----------------------------------------------------------

## Combinaciones de triadas

combinaciones <- gtools::combinations(n = 6,
                                      r = 3,
                                      v = variables,
                                      repeats.allowed = FALSE)

posiciones <- gtools::combinations(n = 6,
                                   r = 3,
                                   repeats.allowed = FALSE)

## Densidad

modelo_densidad <- mixexp::MixModel(frame = datos, 
                                    response = "Densidad", 
                                    mixcomps = variables,
                                    model = 1)

car::vif(modelo_densidad)

graphs_densidad <- 17:20 %>% # Cambiar numeración para las otras triadas
  purrr::map(function(r){
    
    p <- custom_MixturePlot(x = datos[, posiciones[r, 3]], 
                       y = datos[, posiciones[r, 2]], 
                       z = datos[, posiciones[r, 1]],
                       w = modelo_densidad$fitted.values, 
                       x1lab = "",
                       x2lab = "",
                       x3lab = "", 
                       mod = 1, 
                       n.breaks = 6,
                       ck = list(height = 0.70, space = "right"), 
                       corner.labs= combinaciones[r,],
                       cols = TRUE,
                       color.palette = hcl.colors,
                       despts = FALSE,
                       contrs = TRUE)
    
    return(p)
    
  })

do.call(grid.arrange, c(graphs_densidad, ncol = 2))

## Oleosidad

modelo_oleosidad <- mixexp::MixModel(frame = datos, 
                                     response = "Oleosidad", 
                                     mixcomps = variables,
                                     model = 1)

car::vif(modelo_oleosidad)

graphs_oleosidad <- 17:20 %>% # Cambiar numeración para las otras triadas
  purrr::map(function(r){
    
    p <- custom_MixturePlot(x = datos[, posiciones[r, 3]], 
                            y = datos[, posiciones[r, 2]], 
                            z = datos[, posiciones[r, 1]],
                            w = modelo_oleosidad$fitted.values, 
                            x1lab = "",
                            x2lab = "",
                            x3lab = "", 
                            mod = 1, 
                            n.breaks = 6,
                            ck = list(height = 0.70, space = "right"), 
                            corner.labs= combinaciones[r,],
                            cols = TRUE,
                            color.palette = hcl.colors,
                            despts = FALSE,
                            contrs = TRUE)
    
    return(p)
    
  })

do.call(grid.arrange, c(graphs_oleosidad, ncol = 2))

# Función de deseabilidad -------------------------------------------------

## Funciones

densidad <- function(x){
  
  0.6304762*x[1] + 0.9519048*x[2] - 
    0.2195238*x[3] + 0.4519048*x[4] + 
    1.5019048*x[5] + 1.1876190*x[6] 
  
}

oleosidad <- function(x){
  
  2.3573810*x[1] + 5.6930952*x[2] - 
    0.9288095*x[3] + 0.4502381*x[4] + 
    5.2145238*x[5] + 4.5002381*x[6] 
  
}

## Definición de los parámetros para optimización

densidad_des <- desirability::dTarget(low = 0.99,
                                      target = 1.00,
                                      high = 1.01, 
                                      tol = 0)

oleosidad_des <- desirability::dTarget(low = 1.272,
                                       target = 3.00,
                                       high = 4.942)

deseabilidad_compuesta <- desirability::dOverall(densidad_des, 
                                                 oleosidad_des)


## Optimizador

optimizar <- function(x, des_obj){
  
  respuesta_1 <- densidad(x)
  
  respuesta_2 <- oleosidad(x)
  
  prediccion <- predict(des_obj,
                        data.frame(respuesta_1,
                                   respuesta_2))
  
  if(any(abs(x) > 1)){
    
    prediccion <- 0
  }
  
  prediccion
  
}

## Solo correr una vez, para aproximar la región donde se encuentra el 
## óptimo, la siguiente región está acotada a una zona donde tiene
## sentido buscar el óptimo con más detalle

region_experimentacion <- tidyr::expand_grid(x1 = seq(0.05, 0.75, by = 0.20),
                                             x2 = seq(0.05, 0.75, by = 0.20),
                                             x3 = seq(0.05, 0.75, by = 0.20),
                                             x4 = seq(0.05, 0.75, by = 0.20),
                                             x5 = seq(0.05, 0.75, by = 0.20),
                                             x6 = seq(0.05, 0.75, by = 0.20))

for(i in 1:dim(region_experimentacion)[1]) {
  
  tmp <- optim(par = as.vector(region_experimentacion[i,]),
               fn = optimizar,
               des_obj = deseabilidad_compuesta,
               control = list(fnscale = -1))
  if(i == 1)
  {
    best <- tmp
  } else {
    if(tmp$value > best$value && all(tmp$par >= 0.05) && all(tmp$par <= 0.75)) 
      
      best <- tmp
  }
}

round(best$par,2) # combinación de componentes

optimo <- data.frame(densidad(x = best$par),
                     oleosidad(x = best$par))


## Valores de deseabilidad

predict(deseabilidad_compuesta, optimo, all = T) %>% 
  rbind(c(optimo, NA)) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::remove_rownames() %>% 
  setNames(c("Respuesta esperada", "Deseabilidad")) %>% 
  dplyr::mutate(Variable = c("Densidad",
                             "Oleosidad",
                             "Deseabilidad compuesta"),
                Unidades = c("g/ml", "-", "NA")) %>% 
  dplyr::relocate(Variable, .before = 1) %>% 
  dplyr::relocate(Unidades, .before = Deseabilidad)
