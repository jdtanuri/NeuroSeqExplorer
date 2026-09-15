import os
import json
import sqlite3
import time
import re
import warnings
from pathlib import Path
from datetime import datetime

from dotenv import load_dotenv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import BiopythonWarning
from Bio import Entrez

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")

BASES_VALIDAS = set("ATCG")


def validar_sequencia(sequencia: str):
    """Levanta ValueError se a sequência for inválida."""
    if not sequencia:
        raise ValueError("Sequência vazia.")

    sequencia = sequencia.upper()
    invalidos = set(sequencia) - BASES_VALIDAS

    if invalidos:
        raise ValueError(
            f"Sequência contém caractere(s) inválido(s): {', '.join(sorted(invalidos))}"
        )

    return sequencia


def contar_bases(sequencia):
    return {
        "A": sequencia.count("A"),
        "C": sequencia.count("C"),
        "T": sequencia.count("T"),
        "G": sequencia.count("G")
    }


def calcular_gc(sequencia):
    if not sequencia:
        return 0.0

    g = sequencia.count("G")
    c = sequencia.count("C")

    return (g + c) / len(sequencia) * 100


def fita_complementar(sequencia):
    return str(Seq(sequencia).complement())


def transcrever_rna(sequencia):
    return str(Seq(sequencia).transcribe())


def traduzir_proteina(rna):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        return str(Seq(rna).translate(to_stop=True))


def encontrar_orfs_em_fita(sequencia):
    orfs = []

    for frame in range(3):

        trecho = sequencia[frame:]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            rna = str(Seq(trecho).transcribe())

        i = 0
        while i < len(rna) - 2:

            if rna[i:i + 3] != "AUG":
                i += 3
                continue

            rna_orf = rna[i:]

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", BiopythonWarning)
                seq_obj = Seq(rna_orf)
                proteina_completa = str(seq_obj.translate(to_stop=False))
                proteina_ate_stop = str(seq_obj.translate(to_stop=True))

            tem_stop = len(proteina_completa) > len(proteina_ate_stop)

            if proteina_ate_stop:
                orfs.append({
                    "frame": frame + 1,
                    "posicao_inicio": i + frame,
                    "proteina": proteina_ate_stop,
                    "completo": tem_stop
                })

            i += 3

    return orfs


def encontrar_orfs(sequencia):
    orfs = []

    orfs.extend(encontrar_orfs_em_fita(sequencia))

    reversa = str(Seq(sequencia).reverse_complement())

    for orf in encontrar_orfs_em_fita(reversa):
        orf["fita"] = "reversa"
        orfs.append(orf)

    for orf in orfs:
        orf.setdefault("fita", "direta")

    return orfs


def analisar_sequencia(sequencia):

    contagem = contar_bases(sequencia)
    gc = calcular_gc(sequencia)
    complementar = fita_complementar(sequencia)
    rna = transcrever_rna(sequencia)
    proteina = traduzir_proteina(rna)
    orfs = encontrar_orfs(sequencia)

    return {
        "gc": gc,
        "fita_complementar": complementar,
        "rna": rna,
        "proteina": proteina,
        "orfs": orfs,
        **contagem
    }


def ler_fasta(arquivo):
    try:
        registro = SeqIO.read(arquivo.file, "fasta")
    except ValueError as e:
        raise ValueError(f"Arquivo FASTA inválido: {e}")

    return str(registro.seq)


def buscar_gene_ncbi(nome_gene):

    termo_busca = (
        f"{nome_gene}[Gene Name] AND Homo sapiens[Organism] "
        f"AND biomol_mrna[PROP] AND RefSeq[Filter]"
    )

    resultado_busca = Entrez.esearch(
        db="nucleotide",
        term=termo_busca,
        retmax=1
    )
    dados_busca = Entrez.read(resultado_busca)
    resultado_busca.close()

    ids_encontrados = dados_busca["IdList"]

    if not ids_encontrados:
        raise ValueError(f"Nenhuma sequência de mRNA encontrada para o gene '{nome_gene}'.")

    id_escolhido = ids_encontrados[0]

    handle = Entrez.efetch(
        db="nucleotide",
        id=id_escolhido,
        rettype="fasta",
        retmode="text"
    )
    registro = SeqIO.read(handle, "fasta")
    handle.close()

    return str(registro.seq)


def carregar_genes_neurotransmissores():
    caminho = Path(__file__).parent / "genes_neurotransmissores.json"

    with open(caminho, "r", encoding="utf-8") as arquivo:
        return json.load(arquivo)


def buscar_info_gene(nome_gene):
    genes = carregar_genes_neurotransmissores()
    for categoria in genes.values():
        for gene in categoria:
            if gene["gene"].upper() == nome_gene.upper():
                return gene
    return None


def buscar_artigos_recentes(nome_gene, quantidade=5):
    try:
        busca_gene = Entrez.esearch(
            db="gene",
            term=f"{nome_gene}[sym] AND Homo sapiens[orgn]",
            retmax=1
        )
        dados_gene = Entrez.read(busca_gene)
        busca_gene.close()

        ids_gene = dados_gene["IdList"]
        if not ids_gene:
            return []

        gene_id = ids_gene[0]

        time.sleep(0.4)

        link_resultado = Entrez.elink(dbfrom="gene", db="pubmed", id=gene_id)
        link_dados = Entrez.read(link_resultado)
        link_resultado.close()

        if not link_dados[0]["LinkSetDb"]:
            return []

        pmids = [link["Id"] for link in link_dados[0]["LinkSetDb"][0]["Link"]]

        if not pmids:
            return []

        pmids_para_buscar = pmids[:50]

        time.sleep(0.4)

        resultado_resumo = Entrez.esummary(db="pubmed", id=",".join(pmids_para_buscar))
        resumos = Entrez.read(resultado_resumo)
        resultado_resumo.close()

        artigos = []
        for resumo in resumos:
            pmid = str(resumo.get("Id", ""))
            data_publicacao = resumo.get("PubDate", "")

            match_ano = re.search(r"\d{4}", data_publicacao)
            ano = int(match_ano.group()) if match_ano else 0

            artigos.append({
                "titulo": resumo.get("Title", ""),
                "revista": resumo.get("FullJournalName", ""),
                "data_publicacao": data_publicacao,
                "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "ano": ano
            })

        artigos.sort(key=lambda a: a["ano"], reverse=True)

        for artigo in artigos:
            del artigo["ano"]

        return artigos[:quantidade]

    except Exception:
        return []


def inicializar_banco():
    caminho_banco = Path(__file__).parent / "neuroseq.db"

    conexao = sqlite3.connect(caminho_banco)
    cursor = conexao.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS analises (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            gc REAL,
            proteina TEXT,
            quantidade_orfs INTEGER,
            data_hora TEXT
        )
    """)

    conexao.commit()
    conexao.close()


def salvar_analise(gene, gc, proteina, quantidade_orfs):
    caminho_banco = Path(__file__).parent / "neuroseq.db"

    conexao = sqlite3.connect(caminho_banco)
    cursor = conexao.cursor()

    agora = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cursor.execute(
        "INSERT INTO analises (gene, gc, proteina, quantidade_orfs, data_hora) VALUES (?, ?, ?, ?, ?)",
        (gene, gc, proteina, quantidade_orfs, agora)
    )

    conexao.commit()
    conexao.close()


def listar_historico():
    caminho_banco = Path(__file__).parent / "neuroseq.db"

    conexao = sqlite3.connect(caminho_banco)
    cursor = conexao.cursor()

    cursor.execute("SELECT id, gene, gc, proteina, quantidade_orfs, data_hora FROM analises ORDER BY id DESC")
    linhas = cursor.fetchall()

    conexao.close()

    historico = []
    for linha in linhas:
        historico.append({
            "id": linha[0],
            "gene": linha[1],
            "gc": linha[2],
            "proteina": linha[3],
            "quantidade_orfs": linha[4],
            "data_hora": linha[5]
        })

    return historico


def comparar_sequencias(nome1, sequencia1, nome2, sequencia2):
    resultado1 = analisar_sequencia(sequencia1)
    resultado2 = analisar_sequencia(sequencia2)

    return {
        "sequencia_1": {"nome": nome1, **resultado1},
        "sequencia_2": {"nome": nome2, **resultado2},
        "comparacao": {
            "diferenca_gc": resultado1["gc"] - resultado2["gc"],
            "diferenca_tamanho": len(sequencia1) - len(sequencia2),
            "proteinas_identicas": resultado1["proteina"] == resultado2["proteina"],
            "diferenca_quantidade_orfs": len(resultado1["orfs"]) - len(resultado2["orfs"])
        }
    }