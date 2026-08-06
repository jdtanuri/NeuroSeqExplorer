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
    complemento = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    fita = ""

    for letra in sequencia:
        fita += complemento[letra]

    return fita


def transcrever_rna(sequencia):
    return sequencia.replace("T", "U")


def traduzir_proteina(rna):
    codigo_genetico = {
        "UUU": "F", "UUC": "F",
        "UUA": "L", "UUG": "L",

        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",

        "UAU": "Y", "UAC": "Y",
        "UAA": "*", "UAG": "*",

        "UGU": "C", "UGC": "C",
        "UGA": "*",
        "UGG": "W",

        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",

        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",

        "CAU": "H", "CAC": "H",
        "CAA": "Q", "CAG": "Q",

        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

        "AUU": "I", "AUC": "I", "AUA": "I",
        "AUG": "M",

        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",

        "AAU": "N", "AAC": "N",
        "AAA": "K", "AAG": "K",

        "AGU": "S", "AGC": "S",
        "AGA": "R", "AGG": "R",

        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",

        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",

        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",

        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }

    proteina = ""
    traduzindo = False

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]

        if len(codon) < 3:
            break

        if codon == "AUG":
            traduzindo = True

        if not traduzindo:
            continue

        if codon in ["UAA", "UAG", "UGA"]:
            break

        proteina += codigo_genetico[codon]

    return proteina