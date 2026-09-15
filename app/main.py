from fastapi import FastAPI, UploadFile, File, HTTPException
from pydantic import BaseModel, Field
from fastapi.middleware.cors import CORSMiddleware

from app.services import (
    analisar_sequencia,
    ler_fasta,
    validar_sequencia,
    buscar_gene_ncbi,
    carregar_genes_neurotransmissores,
    buscar_info_gene,
    buscar_artigos_recentes,
    inicializar_banco,
    salvar_analise,
    listar_historico,
    comparar_sequencias
)

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

inicializar_banco()


class Sequencia(BaseModel):
    nome: str
    sequencia: str = Field(..., min_length=1)


@app.get("/")
def home():
    return {"message": "Olá, mundo!"}


@app.get("/ola/{nome}")
def ola(nome):
    return {"message": f"Olá, {nome}!"}


@app.post("/analisar")
def analisar(dados: Sequencia):

    try:
        sequencia = validar_sequencia(dados.sequencia)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    resultado = analisar_sequencia(sequencia)

    return {
        "nome": dados.nome,
        **resultado
    }


@app.post("/analisar-fasta")
def analisar_fasta(arquivo: UploadFile = File(...)):

    try:
        sequencia = ler_fasta(arquivo)
        sequencia = validar_sequencia(sequencia)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    return analisar_sequencia(sequencia)


@app.get("/buscar-gene/{nome}")
def buscar_gene(nome: str):

    try:
        sequencia = buscar_gene_ncbi(nome)
        sequencia = validar_sequencia(sequencia)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    resultado = analisar_sequencia(sequencia)
    info_gene = buscar_info_gene(nome)
    artigos = buscar_artigos_recentes(nome)

    salvar_analise(
        gene=nome,
        gc=resultado["gc"],
        proteina=resultado["proteina"],
        quantidade_orfs=len(resultado["orfs"])
    )

    return {
        "gene": nome,
        "info": info_gene,
        "artigos_recentes": artigos,
        **resultado
    }


@app.get("/genes")
def listar_genes():
    return carregar_genes_neurotransmissores()


@app.get("/historico")
def historico():
    return listar_historico()


@app.get("/comparar-genes/{gene1}/{gene2}")
def comparar_genes(gene1: str, gene2: str):

    try:
        sequencia1 = buscar_gene_ncbi(gene1)
        sequencia1 = validar_sequencia(sequencia1)
        sequencia2 = buscar_gene_ncbi(gene2)
        sequencia2 = validar_sequencia(sequencia2)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    resultado = comparar_sequencias(gene1, sequencia1, gene2, sequencia2)

    salvar_analise(
        gene=gene1,
        gc=resultado["sequencia_1"]["gc"],
        proteina=resultado["sequencia_1"]["proteina"],
        quantidade_orfs=len(resultado["sequencia_1"]["orfs"])
    )
    salvar_analise(
        gene=gene2,
        gc=resultado["sequencia_2"]["gc"],
        proteina=resultado["sequencia_2"]["proteina"],
        quantidade_orfs=len(resultado["sequencia_2"]["orfs"])
    )

    return resultado