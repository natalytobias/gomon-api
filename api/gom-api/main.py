from io import StringIO
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from typing import List, Optional
import pandas as pd
from pandas.errors import ParserError

# Cria a instância da aplicação FastAPI
app = FastAPI()

# Endpoint de verificação de status da API
@app.get("/")
async def home():
    return {"message": "Bem-vindo à sua API. Use o endpoint /upload-data/ para enviar seus dados."}

# Endpoint principal para o upload de dados
@app.post("/upload-data/")
async def processar_dados(
    file: UploadFile = File(...),
    k_initial: int = Form(...),
    k_final: int = Form(...),
    case_id: str = Form(...),
    internal_vars: Optional[List[str]] = Form(None)
):
    """
    Recebe um arquivo CSV, variáveis inteiras e uma lista opcional de strings.
    
    Args:
        file (UploadFile): O arquivo CSV.
        k_initial (int): Variável de início.
        k_final (int): Variável de fim.
        case_id (str): Nome da coluna de identificação no CSV.
        internal_vars (List[str]): Lista de variáveis internas do modelo.
    """
    # Validação do tipo de arquivo
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="O arquivo deve ser um CSV.")
        
    try:
        # Acessa o conteúdo do arquivo
        contents = await file.read()
        csv_data = StringIO(contents.decode('utf-8'))
        
        # Carrega o CSV para um DataFrame do Pandas
        df = pd.read_csv(csv_data)

        # --- Validação da coluna 'case_id' ---
        case_id_column = case_id 
        if case_id_column not in df.columns:
            raise HTTPException(status_code=400, detail=f"O CSV não contém a coluna de identificação '{case_id_column}'.")

        # --- Validação das colunas de 'internal_vars' ---
        if internal_vars: 
            missing_vars = [var for var in internal_vars if var not in df.columns]
            print(missing_vars)
            
            if missing_vars:
                # Retorna uma mensagem clara com os nomes das colunas que faltam
                raise HTTPException(status_code=400, detail=f"As seguintes variáveis de 'internal_vars' não foram encontradas no CSV: {', '.join(missing_vars)}")


        # gom.models[[paste0("K", k)]] <- GoMRcpp(
        #     data.object = file,
        #     initial.K = k_initial, final.K = k_final,
        #     gamma.algorithm = "gradient.1992",
        #     initial.gamma = "equal.values",
        #     gamma.fit = TRUE,
        #     lambda.algorithm = "gradient.1992",
        #     initial.lambda = "random",
        #     lambda.fit = TRUE,
        #     case.id = "SubjID",
        #     internal.var = internal_vars
        #     order.K = TRUE,
        #     dec.char = "."
        # ) 

        return {
            "status": "sucesso",
            "message": "Dados recebidos e validados com sucesso!",
            "file_name": file.filename,
            "k_initial": k_initial,
            "k_final": k_final,
            "case_id_column_name": case_id,
            "internal_vars": internal_vars
        }

    except ParserError:
        raise HTTPException(status_code=400, detail="O arquivo CSV está mal formatado e não pôde ser lido.")
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Não foi possível decodificar o arquivo CSV. Verifique a codificação (ex: UTF-8).")
    except Exception as e:
        # Captura e retorna erros inesperados
        raise HTTPException(status_code=500, detail=f"Ocorreu um erro inesperado: {str(e)}")