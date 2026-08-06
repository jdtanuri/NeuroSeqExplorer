from fastapi import FastAPI
from pydantic import BaseModel

from app.services import (
    contar_bases,
    calcular_gc,
    fita_complementar,
    transcrever_rna,
    traduzir_proteina
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

    contagem = contar_bases(sequencia)
    gc = calcular_gc(sequencia)
    complementar = fita_complementar(sequencia)
    rna = transcrever_rna(sequencia)
    proteina = traduzir_proteina(rna)

    return {
        "nome": dados.nome,
        "gc": gc,
        "fita_complementar": complementar,
        "rna": rna,
        "proteina": proteina,
        **contagem
    }