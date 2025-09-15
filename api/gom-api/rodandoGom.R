# Criar dados fictícios compatíveis com o GoM
set.seed(123)

n <- 1000

SubjID <- paste0("ID", 1:n)
Var1 <- factor(sample(1:3, n, replace = TRUE))  # 3 níveis
Var2 <- factor(sample(1:2, n, replace = TRUE))  # 2 níveis
Var3 <- factor(sample(1:4, n, replace = TRUE))  # 4 níveis

teste <- data.frame(SubjID, Var1, Var2, Var3)

##Salva o arquivo
write.csv(teste, "teste.csv", row.names = FALSE)


## Carregando os dados ##
teste <- read.csv("teste.csv", stringsAsFactors = TRUE)

str(teste)
levels(teste$Var1)  # Deve ser: "1" "2" "3"
as.numeric(levels(teste$Var1))  # Deve dar: 1 2 3

#### Rodar o modelo Grade of Membership (GoM) usando Rcpp ####

# Pacotes
#if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
if (!requireNamespace("inline", quietly = TRUE)) install.packages("inline")

#library(readxl)
library(Rcpp)
library(inline)

# Carregar o script GoMRcpp (ajustado para versão atual do R)
source("GoMRcpp.R")


# Executar o GoM para um número específico de perfis
# (exemplo: de 2 a 4 perfis)

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
    case.id = "SubjID",
    internal.var = c("Var1", "Var2", "Var3"),
    order.K = TRUE,
    dec.char = "."
  )
}
