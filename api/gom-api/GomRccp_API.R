#!/usr/bin/env Rscript

# Este script recebe o caminho de um arquivo CSV e parâmetros da linha de comando,
# executa o modelo GoM e salva o resultado em um arquivo de saída.

# Pacotes necessários
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
if (!requireNamespace("inline", quietly = TRUE)) install.packages("inline")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(Rcpp)
library(inline)
library(jsonlite)

# 1. Função para analisar os argumentos de linha de comando
# Esta função converte os argumentos "--chave valor" em uma lista nomeada.
parse_args <- function(args) {
  params <- list()
  for (i in seq_along(args)) {
    if (startsWith(args[i], "--")) {
      key <- substr(args[i], 3, nchar(args[i]))
      if (i + 1 <= length(args) && !startsWith(args[i+1], "--")) {
        params[[key]] <- args[i+1]
      } else {
        params[[key]] <- TRUE
      }
    }
  }
  return(params)
}

# 2. Pega e analisa os argumentos da linha de comando
args <- commandArgs(trailingOnly = TRUE)
params <- parse_args(args)

# Verificação básica dos parâmetros
required_params <- c("file-path", "k-initial", "k-final", "case-id")
if (!all(required_params %in% names(params))) {
  stop("Parâmetros obrigatórios ausentes: --file-path, --k-initial, --k-final, --case-id")
}

# Extrai os parâmetros
file_path <- params[["file-path"]]
k_initial <- as.integer(params[["k-initial"]])
k_final <- as.integer(params[["k-final"]])
case_id_column <- params[["case-id"]]
# Converte a string de variáveis internas de volta para um vetor, se existir
internal_vars_str <- params[["internal-vars"]]
internal_vars <- if (!is.null(internal_vars_str) && nchar(internal_vars_str) > 0) unlist(strsplit(internal_vars_str, ",")) else NULL
output_path <- params[["output-path"]]

# 3. Carrega o script GoMRcpp.R que contém a função principal
if (file.exists("GoMRcpp.R")) {
  source("GoMRcpp.R")
} else {
  stop("O arquivo GoMRcpp.R não foi encontrado no diretório do projeto.")
}

# 4. Lê os dados do CSV fornecido
if (!file.exists(file_path)) {
  stop(paste("Arquivo não encontrado:", file_path))
}
data_object <- read.csv(file_path, stringsAsFactors = TRUE)

# 5. Validações adicionais (opcional, mas recomendado)
if (case_id_column %in% names(data_object)) {
  data_object[[case_id_column]] <- as.factor(data_object[[case_id_column]])
} else {
  stop(paste("A coluna de identificação '", case_id_column, "' não foi encontrada no CSV.", sep=""))
}

if (!is.null(internal_vars)) {
  missing_vars <- setdiff(internal_vars, names(data_object))
  if (length(missing_vars) > 0) {
    stop(paste("As seguintes variáveis internas não foram encontradas no CSV:", paste(missing_vars, collapse=", ")))
  }
}

# 6. Executa o modelo GoM para cada valor de K no intervalo
# Aqui você pode adaptar o loop conforme a sua necessidade.
# O exemplo a seguir roda o modelo apenas para um K específico.
# Para rodar para um intervalo, você usaria um loop.
cat(paste("Executando modelo com K de", k_initial, "a", k_final, "perfis...\n"))

# Adaptado para o seu loop original
gom.models <- list()
for (k in k_initial:k_final) {
  cat(paste("Executando modelo com", k, "perfis...\n"))
  
  gom.models[[paste0("K", k)]] <- GoMRcpp(
    data.object = data_object,
    initial.K = k, final.K = k,
    gamma.algorithm = "gradient.1992",
    initial.gamma = "equal.values",
    gamma.fit = TRUE,
    lambda.algorithm = "gradient.1992",
    initial.lambda = "random",
    lambda.fit = TRUE,
    case.id = case_id_column,
    internal.var = internal_vars,
    order.K = TRUE,
    dec.char = "."
  )
}

# 7. Salva o resultado em um arquivo JSON para que o Python possa lê-lo
# Isso é crucial para que a sua API possa retornar os resultados.
json_output <- jsonlite::toJSON(gom.models, pretty = TRUE)
write(json_output, file = output_path)

# Encerra o script com sucesso
quit(save = "no", status = 0)