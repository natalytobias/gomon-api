import subprocess
import os
import tempfile
import json
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
    """
    Retorna uma mensagem de boas-vindas para indicar que a API está funcionando.
    """
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
    
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="O arquivo deve ser um CSV.")
        
    try:
        # Acessa o conteúdo do arquivo para validação inicial (sem salvar no disco)
        contents = await file.read()
        csv_data = StringIO(contents.decode('utf-8'))
        df = pd.read_csv(csv_data)

        # Validação das colunas de entrada no DataFrame
        if case_id not in df.columns:
            raise HTTPException(status_code=400, detail=f"O CSV não contém a coluna de identificação '{case_id}'.")
        
        if internal_vars:
            missing_vars = [var for var in internal_vars if var not in df.columns]
            if missing_vars:
                raise HTTPException(status_code=400, detail=f"As seguintes variáveis de 'internal_vars' não foram encontradas no CSV: {', '.join(missing_vars)}")

        # --- Chamada ao script R com subprocess ---
        # Cria um diretório temporário para salvar o CSV e o arquivo de saída
        with tempfile.TemporaryDirectory() as temp_dir:
            csv_path = os.path.join(temp_dir, file.filename)
            
            # Reposiciona o ponteiro do arquivo para o início para salvá-lo
            await file.seek(0)
            with open(csv_path, "wb") as f:
                f.write(await file.read())

            # Define o caminho do arquivo de saída JSON
            output_file_path = os.path.join(temp_dir, "model_output.json")
            
            # Converte a lista de variáveis internas para uma string separada por vírgulas
            internal_vars_str = ",".join(internal_vars) if internal_vars else ""

            # Monta a lista de argumentos para a chamada de subprocesso
            cmd_args = [
                "Rscript",
                "GomRccp_API.R",
                "--file-path", csv_path,
                "--k-initial", str(k_initial),
                "--k-final", str(k_final),
                "--case-id", case_id,
                "--output-path", output_file_path
            ]
            
            # Adiciona as variáveis internas apenas se elas existirem
            if internal_vars_str:
                cmd_args.extend(["--internal-vars", internal_vars_str])

            # Executa o script R
            result = subprocess.run(cmd_args, capture_output=True, text=True)

            if result.returncode != 0:
                # Se o script R falhar, retorna o erro gerado pelo R
                raise HTTPException(status_code=500, detail=f"Erro no script R: {result.stderr}")

            # 7. Lê e retorna o resultado do arquivo de saída JSON
            if not os.path.exists(output_file_path):
                raise HTTPException(status_code=500, detail="O script R não gerou o arquivo de saída esperado.")
            
            with open(output_file_path, "r") as f:
                r_output = json.load(f)
            
            return {
                "status": "sucesso",
                "message": "Dados processados pelo script R com sucesso!",
                "r_output": r_output,
                "file_name": file.filename
            }

    except ParserError:
        raise HTTPException(status_code=400, detail="O arquivo CSV está mal formatado e não pôde ser lido.")
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Não foi possível decodificar o arquivo CSV. Verifique a codificação (ex: UTF-8).")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Ocorreu um erro inesperado: {str(e)}")