from Bio.Seq import Seq
from Bio import SeqIO


def contar_bases(sequencia):
    return {
        "A": sequencia.count("A"),
        "C": sequencia.count("C"),
        "T": sequencia.count("T"),
        "G": sequencia.count("G")
    }


def calcular_gc(sequencia):
    g = sequencia.count("G")
    c = sequencia.count("C")

    return (g + c) / len(sequencia) * 100


def fita_complementar(sequencia):
    return str(Seq(sequencia).complement())


def transcrever_rna(sequencia):
    return str(Seq(sequencia).transcribe())


def traduzir_proteina(rna):
    return str(Seq(rna).translate(to_stop=True))


def encontrar_orfs_em_fita(sequencia):
    orfs = []

    for frame in range(3):

        trecho = sequencia[frame:]

        rna = str(Seq(trecho).transcribe())

        inicio = -1

        for i in range(0, len(rna), 3):
            if rna[i:i+3] == "AUG":
                inicio = i
                break

        if inicio == -1:
            continue

        rna_orf = rna[inicio:]

        proteina = str(Seq(rna_orf).translate(to_stop=True))

        if proteina:
            orfs.append({
                "frame": frame + 1,
                "proteina": proteina
            })

    return orfs


def encontrar_orfs(sequencia):
    orfs = []

    orfs.extend(encontrar_orfs_em_fita(sequencia))

    reversa = str(Seq(sequencia).reverse_complement())

    orfs.extend(encontrar_orfs_em_fita(reversa))

    return orfs


def analisar_sequencia(sequencia):

    contagem = contar_bases(sequencia)
    gc = calcular_gc(sequencia)
    complementar = fita_complementar(sequencia)
    rna = transcrever_rna(sequencia)
    proteina = traduzir_proteina(rna)

    return {
        "gc": gc,
        "fita_complementar": complementar,
        "rna": rna,
        "proteina": proteina,
        **contagem
    }


def ler_fasta(arquivo):
    registro = SeqIO.read(arquivo.file, "fasta")

    return str(registro.seq)