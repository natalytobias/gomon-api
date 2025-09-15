# Pacotes
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
if (!requireNamespace("inline", quietly = TRUE)) install.packages("inline")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr") 

library(Rcpp)
library(inline)
library(readr)

# Carregar o script GoMRcpp (ajustado para versão atual do R)
source("GoMRcpp.R")

# Carregar os dados do arquivo "clicamama.csv"
teste <- read_csv("clicamama.csv")

# 1. Defina o ID de cada caso
# A coluna "genename" é o identificador único e não deve ser usada na análise.
id_paciente <- "genename"

# 2. Selecione TODAS as outras variáveis internas para a análise.
# Usamos a função names() para pegar todas as colunas do data frame e removemos "genename".
variaveis_gom <- setdiff(names(teste), id_paciente)

#### Rodar o modelo Grade of Membership (GoM) usando Rcpp ####

gom.models <- list()  # para armazenar os modelos

for (k in 2:4) {
  cat(paste("Executando modelo com", k, "perfis...\n"))
  
  gom.models[[paste0("K", k)]] <- GoMRcpp(
    data.object = teste,
    initial.K = k, final.K = k,
    gamma.algorithm = "gradient.1992",
    initial.gamma = "equal.values",
    gamma.fit = TRUE,
    lambda.algorithm = "gradient.1992",
    initial.lambda = "random",
    lambda.fit = TRUE,
    case.id = id_paciente,
    internal.var = variaveis_gom, # Usando todas as variáveis restantes
    order.K = TRUE,
    dec.char = "."
  )
}