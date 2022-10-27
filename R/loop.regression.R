#' @title Performs combined bivariate and multivariate linear regressions.
#' @description The loop.regression function is used to perform loop fits of linear models, including multivariate .
#' @param data the data parameter is a data.frame.
#' @param varY is a character vector containing the names of the dependent variables.
#' @param varX is a character vector containing the names of the independent variables.
#' @param varXadj is a character vector with the names of the adjustment variables.
#' @param varXconf is a character vector with the names of the confounders.
#' @return Se obtiene un data.frame con los siguientes datos: varY, varX, estimate, std.error, conf.level, conf.low, conf.high, statistic, df.error, p.value, N, call (variables usadas), model (unadjusted, adjusted), type (lineal, logistic)
#' @examples
#' table <- loop.regression(data = data, varY = c("varY_1", "varY_2"), varX = c("varX_1", "varX_2"), varXadj = "creatinine",
#'    varXconf = c("MaternalAge", "GestationalAge"),  creatinine=TRUE, multivariate=TRUE)
#' table <- loop.regression(iris, varY = "Sepal.Length", varX = c("Species","Sepal.Width"), varXadj = "Petal.Length",
#'    varXconf = "Petal.Width", creatinine = TRUE, multivariate = TRUE)


#Description
## La función loop.regression es usada para realizar bucles de ajustes de modelos lineales, incluyendo multivariados

# Usage
## loop.regression(data, varY, varX, varXadj, varXconf, creatinine=FALSE, multivariate=FALSE)
## data: es un data frame que contiene las varibles del modelo
## varY: es un vector de caracteres con los nombres de las variables depedientes. varY: Si la variable Y es numérica hace regresión
##    lineal y si es factorial hace regresión logística y calcula el OR y sus intervalos de confianza.
## varX: es un vector de caracteres con los nombres de las variables independientes. Si la variable X es factorial, reconoce cuántas
##    categorías tiene para exportar todas a la tabla menos la de referencia
## varXadj: es un vector de caracteres con los nombres de las variables de ajuste. Podemos indicar que ajuste por creatinina
##    (creatinine = TRUE) y le indicamos varXadj.
## varXconf: es un vector de caracteres con los nombres de las variables confusoras. Podemos indecarle que queremos que haga, además
##    de bivariado que hace por defecto, que también haga el multivariado (multivariate = TRUE) y le indicamos qué covariables usar (varXconf).

# Arguments
## creatinine: logical. Por defecto es FALSE y no usaría la variable varXadj, pero si es TRUE añade al modelo la variable varXadj
## multivariate: logical. Por defecto es FALSE y ajusta el modelo lineal bivariado, pero si es TRUE ajusta el modelo lineal
##    multivariado añadiendo como confusores el vector varXconf

#Details
## Obtenemos una tabla con los siguientes datos: varY, varX, estimate, std.error, conf.level, conf.low, conf.high, statistic,
##    df.error, p.value, N, call (variables usadas), model (unadjusted, adjusted), type (lineal, logistic)


# Examples:
## require(broomExtra)
## table <- loop.regression(data = data, varY = c("varY_1", "varY_2"), varX = c("varX_1", "varX_2"), varXadj = "creatinine",
##    varXconf = c("MaternalAge", "GestationalAge"),  creatinine=TRUE, multivariate=TRUE)
## table <- loop.regression(iris, varY = "Sepal.Length", varX = c("Species","Sepal.Width"), varXadj = "Petal.Length",
##    varXconf = "Petal.Width", creatinine = TRUE, multivariate = TRUE)

loop.regression <- function(data, varY, varX, varXadj, varXconf, creatinine=FALSE, multivariate=FALSE) {

require(broomExtra)

if(creatinine==TRUE && multivariate==TRUE){

  #bivariate
  for(j in 1:length(varY)){
    for(i in 1:length(varX)){

      if(is.numeric(data[varY[j]][,1])==TRUE) {

        modelo <- lm(reformulate(c(varX[i], varXadj),varY[j]), data = data, na.action = na.exclude)
        call <- deparse(formula(modelo), width.cutoff = 500L)
        N <- nrow(modelo$model)
        row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
        model <- "unadjusted"
        type <- "lineal"

        if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
          table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
        }
        else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
          table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
        }

      } else{

        modelo <- glm(reformulate(c(varX[i], varXadj),varY[j]), family="binomial", data = data, na.action = na.exclude)
        call <- deparse(formula(modelo), width.cutoff = 500L)
        N <- nrow(modelo$model)
        row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
        model <- "unadjusted"
        type <- "logistic"

        if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
          table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
        }
        else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
          table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
        }
      }
    }
  }

  table_model_1 <- table_model

  #multivariate
  for(j in 1:length(varY)){
    for(i in 1:length(varX)){

      if(is.numeric(data[varY[j]][,1])==TRUE) {

        modelo <- lm(reformulate(c(varX[i], varXadj, varXconf),varY[j]), data = data, na.action = na.exclude)
        call <- deparse(formula(modelo), width.cutoff = 500L)
        N <- nrow(modelo$model)
        row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
        model <- "adjusted"
        type <- "lineal"

        if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
          table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
        }
        else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
          table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
        }

      } else{

        modelo <- glm(reformulate(c(varX[i], varXadj, varXconf),varY[j]), family="binomial", data = data, na.action = na.exclude)
        call <- deparse(formula(modelo), width.cutoff = 500L)
        N <- nrow(modelo$model)
        row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
        model <- "adjusted"
        type <- "logistic"

        if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
          table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
        }
        else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
          table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
        }
      }
    }
  }

 table_model_1 <- rbind(table_model_1, table_model)

} else {

if(creatinine==FALSE && multivariate==TRUE){

    #bivariate
    for(j in 1:length(varY)){
      for(i in 1:length(varX)){

        if(is.numeric(data[varY[j]][,1])==TRUE) {

          modelo <- lm(reformulate(c(varX[i]),varY[j]), data = data, na.action = na.exclude)
          call <- deparse(formula(modelo), width.cutoff = 500L)
          N <- nrow(modelo$model)
          row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
          model <- "unadjusted"
          type <- "lineal"

          if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
            table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
          }
          else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
            table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
          }

        } else{

          modelo <- glm(reformulate(c(varX[i]),varY[j]), family="binomial", data = data, na.action = na.exclude)
          call <- deparse(formula(modelo), width.cutoff = 500L)
          N <- nrow(modelo$model)
          row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
          model <- "unadjusted"
          type <- "logistic"

          if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
            table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
          }
          else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
            table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
          }
        }
      }
    }

    table_model_1 <- table_model

    #multivariate
    for(j in 1:length(varY)){
      for(i in 1:length(varX)){

        if(is.numeric(data[varY[j]][,1])==TRUE) {

          modelo <- lm(reformulate(c(varX[i], varXconf),varY[j]), data = data, na.action = na.exclude)
          call <- deparse(formula(modelo), width.cutoff = 500L)
          N <- nrow(modelo$model)
          row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
          model <- "adjusted"
          type <- "lineal"

          if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
            table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
          }
          else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
            table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
          }

        } else{

          modelo <- glm(reformulate(c(varX[i], varXconf),varY[j]), family="binomial", data = data, na.action = na.exclude)
          call <- deparse(formula(modelo), width.cutoff = 500L)
          N <- nrow(modelo$model)
          row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
          model <- "adjusted"
          type <- "logistic"

          if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
            table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
          }
          else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
            table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
          }
        }
      }
    }

  table_model_1 <- rbind(table_model_1, table_model)

  } else {

if(creatinine==TRUE && multivariate==FALSE){

      #bivariate
      for(j in 1:length(varY)){
        for(i in 1:length(varX)){

          if(is.numeric(data[varY[j]][,1])==TRUE) {

            modelo <- lm(reformulate(c(varX[i], varXadj),varY[j]), data = data, na.action = na.exclude)
            call <- deparse(formula(modelo), width.cutoff = 500L)
            N <- nrow(modelo$model)
            row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
            model <- "unadjusted"
            type <- "lineal"

            if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
              table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
            }
            else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
              table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
            }

          } else{

            modelo <- glm(reformulate(c(varX[i], varXadj),varY[j]), family="binomial", data = data, na.action = na.exclude)
            call <- deparse(formula(modelo), width.cutoff = 500L)
            N <- nrow(modelo$model)
            row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
            model <- "unadjusted"
            type <- "logistic"

            if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
              table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
            }
            else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
              table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
            }
          }
        }
      }

      table_model_1 <- table_model

  } else {

if(creatinine==FALSE && multivariate==FALSE){

      #bivariate
      for(j in 1:length(varY)){
        for(i in 1:length(varX)){

          if(is.numeric(data[varY[j]][,1])==TRUE) {

            modelo <- lm(reformulate(c(varX[i]),varY[j]), data = data, na.action = na.exclude)
            call <- deparse(formula(modelo), width.cutoff = 500L)
            N <- nrow(modelo$model)
            row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
            model <- "unadjusted"
            type <- "lineal"

            if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
              table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
            }
            else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
              table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
            }

          } else{

            modelo <- glm(reformulate(c(varX[i]),varY[j]), family="binomial", data = data, na.action = na.exclude)
            call <- deparse(formula(modelo), width.cutoff = 500L)
            N <- nrow(modelo$model)
            row <- if(is.numeric(data[varX[i]][,1])==TRUE | length(levels(data[varX[i]][,1]))<3) {2} else{length(levels(data[varX[i]][,1]))}
            model <- "unadjusted"
            type <- "logistic"

            if(i == 1 && j == 1){    # En el primer modelo creamos la tabla, con la primera fila
              table_model <- data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type)
            }
            else{     # En los siguientes modelos la tabla ya existe, y solo hay que a?dirle filas por debajo
              table_model <- rbind(table_model, data.frame(varY[j], tidy_parameters(modelo)[2:row,], N, call, model, type))
            }
          }
        }
      }
  table_model_1 <- table_model
    }

    }
  }
}

  ######

  table_model_1$estimate <- ifelse(table_model_1$type == "logistic", exp(table_model_1$estimate), table_model_1$estimate)
  table_model_1$conf.low <- ifelse(table_model_1$type == "logistic", exp(table_model_1$conf.low), table_model_1$conf.low)
  table_model_1$conf.high <- ifelse(table_model_1$type == "logistic", exp(table_model_1$conf.high), table_model_1$conf.high)

  names(table_model_1)[names(table_model_1) == "varY.j."] <- "varY"
  names(table_model_1)[names(table_model_1) == "term"] <- "varX"

  return(table_model_1)
}
