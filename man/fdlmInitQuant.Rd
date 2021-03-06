% Generated by roxygen2 (4.0.1): do not edit by hand
\name{fdlmInitQuant}
\alias{fdlmInitQuant}
\title{Quantidades iniciais para Gibbs F-DLM}
\usage{
fdlmInitQuant(model, modelOut = TRUE)
}
\arguments{
\item{model}{modelo F-DLM com valores da matriz de dados, de exogenas e hiperparametros da priori.}
}
\value{
Lista com modelo e valores iniciais para a cadeia do Gibbs.
}
\description{
Calculo de quantidades inicias para evitar repeticao de operacoes identicas a cada iteracao do algoritmo de Gibbs..
}

