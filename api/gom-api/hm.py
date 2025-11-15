from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from typing import Optional, List, Dict, Any
import csv
import json
from io import StringIO
import pandas as pd

app = FastAPI()

@app.post("/convertendo-csv")
async def csv_to_json(
    num_k: int = Form(..., description="Número de clusters K"),
    internal_vars_string: Optional[str] = Form(None, description="Variáveis internas"),
    file: UploadFile = File(..., description="Arquivo CSV no formato LMFR")
):
    """
    Converte arquivos CSV no formato LMFR para JSON estruturado
    """
    try:
        # Validar se o arquivo é CSV
        if not file.filename.lower().endswith('.csv'):
            raise HTTPException(
                status_code=400, 
                detail="O arquivo deve ser um CSV (.csv)"
            )
        
        # Ler o conteúdo do arquivo
        contents = await file.read()
        csv_content = contents.decode('utf-8')
        
        # Converter CSV para JSON estruturado
        resultado = processar_csv_lmfr(csv_content, num_k, internal_vars_string)
        
        return resultado
        
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"Erro na conversão: {str(e)}"
        )

def processar_csv_lmfr(csv_content: str, num_k: int, internal_vars_string: Optional[str]) -> Dict[str, Any]:
    """
    Processa arquivos CSV no formato LMFR e converte para JSON estruturado
    """
    dados = []
    csv_file = StringIO(csv_content)
    
    # Ler CSV
    reader = csv.DictReader(csv_file)
    
    for linha in reader:
        # Processar cada linha
        linha_processada = processar_linha_lmfr(linha)
        dados.append(linha_processada)
    
    # Estruturar os dados por variável
    dados_estruturados = estruturar_por_variavel(dados)
    
    # Calcular estatísticas
    estatisticas = calcular_estatisticas(dados)
    
    return {
        "num_k": num_k,
        "internal_vars_string": internal_vars_string,
        "estatisticas": estatisticas,
        "variaveis": dados_estruturados,
        "dados_brutos": dados,
        "total_registros": len(dados),
        "variaveis_unicas": list(dados_estruturados.keys())
    }

def processar_linha_lmfr(linha: Dict[str, str]) -> Dict[str, Any]:
    """
    Processa uma linha do CSV LMFR e converte tipos de dados
    """
    linha_processada = {}
    
    for chave, valor in linha.items():
        # Converter para tipos apropriados
        if chave in ['n', 'perc']:
            # n e perc são numéricos
            linha_processada[chave] = float(valor) if valor else 0.0
        elif chave.startswith('k') and not chave.endswith('_perc_lj'):
            # Valores k são numéricos
            linha_processada[chave] = float(valor) if valor else 0.0
        elif chave.endswith('_perc_lj'):
            # Percentuais LJ são numéricos
            linha_processada[chave] = float(valor) if valor else 0.0
        else:
            # Manter como string
            linha_processada[chave] = valor
    
    return linha_processada

def estruturar_por_variavel(dados: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Estrutura os dados agrupando por variável
    """
    variaveis = {}
    
    for item in dados:
        var_name = item['Variable']
        level = item['Level']
        
        if var_name not in variaveis:
            variaveis[var_name] = {
                "nome": var_name,
                "levels": [],
                "total_n": 0,
                "distribuicao": []
            }
        
        # Adicionar level
        level_data = {
            "level": level,
            "n": item['n'],
            "perc": item['perc'],
            "valores_k": {},
            "percentuais_lj": {}
        }
        
        # Extrair valores k
        for chave, valor in item.items():
            if chave.startswith('k') and not chave.endswith('_perc_lj') and chave != 'k':
                level_data["valores_k"][chave] = valor
            elif chave.endswith('_perc_lj'):
                level_data["percentuais_lj"][chave] = valor
        
        variaveis[var_name]["levels"].append(level_data)
        variaveis[var_name]["total_n"] += item['n']
        variaveis[var_name]["distribuicao"].append({
            "level": level,
            "n": item['n'],
            "perc": item['perc']
        })
    
    return variaveis

def calcular_estatisticas(dados: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calcula estatísticas gerais dos dados
    """
    total_n = sum(item['n'] for item in dados)
    variaveis_unicas = list(set(item['Variable'] for item in dados))
    levels_por_variavel = {}
    
    for item in dados:
        var_name = item['Variable']
        if var_name not in levels_por_variavel:
            levels_por_variavel[var_name] = []
        if item['Level'] not in levels_por_variavel[var_name]:
            levels_por_variavel[var_name].append(item['Level'])
    
    # Encontrar colunas k disponíveis
    colunas_k = set()
    colunas_perc_lj = set()
    
    for item in dados:
        for chave in item.keys():
            if chave.startswith('k') and not chave.endswith('_perc_lj'):
                colunas_k.add(chave)
            elif chave.endswith('_perc_lj'):
                colunas_perc_lj.add(chave)
    
    return {
        "total_observacoes": total_n,
        "numero_variaveis": len(variaveis_unicas),
        "variaveis": variaveis_unicas,
        "levels_por_variavel": levels_por_variavel,
        "colunas_k_disponiveis": sorted(list(colunas_k)),
        "colunas_perc_lj_disponiveis": sorted(list(colunas_perc_lj)),
        "total_registros_csv": len(dados)
    }

# Versão alternativa usando pandas para processamento mais eficiente
@app.post("/convertendo-csv-pandas")
async def csv_to_json_pandas(
    num_k: int = Form(...),
    internal_vars_string: Optional[str] = Form(None),
    file: UploadFile = File(...)
):
    """
    Converte arquivos LMFR para JSON usando pandas
    """
    try:
        # Ler arquivo com pandas
        contents = await file.read()
        csv_content = contents.decode('utf-8')
        csv_file = StringIO(csv_content)
        
        df = pd.read_csv(csv_file)
        
        # Converter para dicionário
        dados = df.to_dict('records')
        
        # Processar dados
        resultado = processar_com_pandas(df, num_k, internal_vars_string)
        
        return resultado
        
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"Erro na conversão: {str(e)}"
        )

def processar_com_pandas(df: pd.DataFrame, num_k: int, internal_vars_string: Optional[str]) -> Dict[str, Any]:
    """
    Processa os dados usando pandas
    """
    # Estatísticas básicas
    estatisticas = {
        "total_observacoes": df['n'].sum(),
        "numero_variaveis": df['Variable'].nunique(),
        "variaveis_unicas": df['Variable'].unique().tolist(),
        "total_registros": len(df)
    }
    
    # Estruturar por variável
    variaveis = {}
    
    for var_name in df['Variable'].unique():
        var_data = df[df['Variable'] == var_name]
        
        levels_data = []
        for _, row in var_data.iterrows():
            level_info = {
                "level": row['Level'],
                "n": float(row['n']),
                "perc": float(row['perc']),
                "valores_k": {},
                "percentuais_lj": {}
            }
            
            # Extrair colunas k
            k_cols = [col for col in df.columns if col.startswith('k') and not col.endswith('_perc_lj')]
            for col in k_cols:
                if col in row:
                    level_info["valores_k"][col] = float(row[col]) if pd.notna(row[col]) else 0.0
            
            # Extrair percentuais LJ
            perc_cols = [col for col in df.columns if col.endswith('_perc_lj')]
            for col in perc_cols:
                if col in row:
                    level_info["percentuais_lj"][col] = float(row[col]) if pd.notna(row[col]) else 0.0
            
            levels_data.append(level_info)
        
        variaveis[var_name] = {
            "nome": var_name,
            "total_n": var_data['n'].sum(),
            "levels": levels_data,
            "distribuicao": var_data[['Level', 'n', 'perc']].to_dict('records')
        }
    
    return {
        "num_k": num_k,
        "internal_vars_string": internal_vars_string,
        "estatisticas": estatisticas,
        "variaveis": variaveis,
        "metadados": {
            "formato": "LMFR",
            "colunas_originais": df.columns.tolist(),
            "colunas_k": [col for col in df.columns if col.startswith('k') and not col.endswith('_perc_lj')],
            "colunas_perc_lj": [col for col in df.columns if col.endswith('_perc_lj')]
        }
    }