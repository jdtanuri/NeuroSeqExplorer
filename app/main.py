from fastapi import FastAPI, UploadFile, File
from pydantic import BaseModel

from app.services import (
    analisar_sequencia,
    ler_fasta
)

app = FastAPI()


class Sequencia(BaseModel):
    nome: str
    sequencia: str


@app.get("/")
def home():
    return {"message": "Olá, mundo!"}


@app.get("/ola/{nome}")
def ola(nome):
    return {"message": f"Olá, {nome}!"}


@app.post("/analisar")
def analisar(dados: Sequencia):

    sequencia = dados.sequencia.upper()

    for letra in sequencia:
        if letra not in "ATCG":
            return {
                "erro": "Sequência contém caractere inválido."
            }

    resultado = analisar_sequencia(sequencia)

    return {
        "nome": dados.nome,
        **resultado
    }


@app.post("/analisar-fasta")
def analisar_fasta(arquivo: UploadFile = File(...)):

    sequencia = ler_fasta(arquivo)

    return analisar_sequencia(sequencia)