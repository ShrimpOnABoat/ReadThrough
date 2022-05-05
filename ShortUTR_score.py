#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:54:34 2019

@author: tony
"""

'''
Principe : 
    1- Pour chaque gène, exon par exon, récupérer la séquence nucléique.
    2- Vérifier la présence d'un codon STOP dans le cadre de lecture.
    3- S'il n'y en a pas, ajouter la séquence suivante et chercher à nouveau un codon STOP in-frame
    4- Lorsqu'un codon STOP est trouvé (ou la fin de la séquence), mémoriser les positions des différents segments,  ainsi que le codon en question
    5- Calculer les sommes et moyennes pour chaque gène.
    6- Générer un fichier contenant le nom du gène, les positions des UTR 3' jusqu'au codon STOP, et les valeurs
'''

from Bio import SeqIO

def find_STOP_codon(sequence):
    # si un codon STOP est présent dans les 30 premiers nucléotides, on l'ignore
    for i in range (30,len(sequence),3):
        if sequence[i:i+3] in ('TAA','TAG','TGA'):
            return i+1, str(sequence[i:i+3])
    return 0, ''

############################################################3
# Ouverture et balayage du fichier UTR3 au format bed

old_gene=''
old_chrom=''
boucles=100
skip_gene=False

fi=open('UTR3_RB_Norm.bed') # Les fichiers UTR3_Shrt_Norm et Sen sont identiques, en fin de compte
fo=open('UTR3_Shrt_codon','w+')

for ligne in fi:
    chaine = (ligne.rstrip('\n')).split('\t')
    # récupération des informations
    chrom = chaine[0]
    start = int(chaine[1])
    end = int(chaine[2])
    gene = chaine[3]
    direction = chaine[5]

    # ouvrir le bon fichier fasta
    if(old_chrom!=chrom):
        file_name = "../../Genome_hg19/Homo_sapiens.GRCh37.dna.chromosome." + chrom[3:] + ".fa"
        chrom_seq=SeqIO.read(file_name, "fasta")
        print("New chromosome : ", chrom)

    # test si nouveau gene
    if(old_gene!=gene):
        skip_gene=False
        # Lecture de la séquence
        # ATTENTION AU SENS !!!!
        if (direction=="-"):
            UTR3=chrom_seq.seq[start:end].reverse_complement()
        else:
            UTR3=chrom_seq.seq[start:end]
        # la séquence contient-elle un codon STOP ?
        stop_codon, codon = find_STOP_codon(UTR3)
        if(stop_codon):
            # codon STOP trouvé, donc on mémorise les positions de début et fin pour ce gène, et on passe au gène suivant, selon la direction
            if(direction=='+'):
                fo.write('\t'.join([chrom, str(start), str(start+stop_codon), gene, direction, codon])+'\n')
            else:
                fo.write('\t'.join([chrom, str(end-stop_codon), str(end), gene, direction, codon])+'\n')
            skip_gene = True
        else:
            # Pas de codon STOP dans cet exon, donc on écrit l'exon et on passe au suivant
            fo.write('\t'.join([chrom, str(start), str(end), gene, direction, codon])+'\n')
            
    elif(not skip_gene):
        # On regarde les exons suivants, car les précédents n'avaient pas de codon STOP
        # Lecture de la séquence et concaténation à la précédente
        if (direction=="-"):
            UTR3+=chrom_seq.seq[start:end].reverse_complement()
        else:
            UTR3+=chrom_seq.seq[start:end]
        # la séquence contient-elle un codon STOP ?
        stop_codon, codon = find_STOP_codon(UTR3)
        if(stop_codon):
            # codon STOP trouvé, donc on mémorise les positions de début et fin pour ce gène, et on passe au gène suivant, selon la direction
            if(direction=='+'):
                fo.write('\t'.join([chrom, str(start), str(start+stop_codon), gene, direction, codon])+'\n')
            else:
                fo.write('\t'.join([chrom, str(end-stop_codon), str(end), gene, direction, codon])+'\n')
            skip_gene = True
        else:
            # Pas de codon STOP dans cet exon, donc on écrit l'exon et on passe au suivant
            fo.write('\t'.join([chrom, str(start), str(end), gene, direction, codon])+'\n')
        

    # mémorisation des dernières valeurs (nécessaire ?)
    old_chrom = chrom
    old_start = start
    old_end = end
    old_gene = gene

#    boucles-=1
#    if (boucles<1):
#        break

fi.close()
fo.close()
