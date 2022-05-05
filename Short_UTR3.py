# Création d'un fichier contenant les résultats du Ribo-Seq pour les UTR seulement, dans le génome hg19

'''
Principe : création d'un tableau contenant les extrémités des UTR, 5' et 3', le chromosome,
le brin, et finalement la moyenne et l'écart-type des valeurs Ribo-Seq
'''
from statistics import mean
import sys


# ouverture des fichiers
NomFichRiboSeq = 'GSM1047589_C_Tx_alignment_reads_coverage_per_position_bothStrands_SF_0.93_.wig'
NomFichUTR5 = 'UTR5'
NomFichUTR3 = 'UTR3'
NomFichKG = 'knownGene.txt'
NomFichXRef = 'kgXref.txt'

'''
# Table XRef : (nom UCSC;nom usuel)
lstXRef = []
lstRef = []
FichXRef = open(NomFichXRef)
for ligne in FichXRef:
    lstXRef.append([ligne.split('\t')[0],ligne.split('\t')[4]])
    lstRef.append(ligne.split('\t')[0])
FichXRef.close()

# Table UTR3 : (chromosome, debut, fin, nom usuel, sens)
lstUTR3 = []
FichUTR3 = open(NomFichUTR3)
for ligne in FichUTR3:
    lstUTR3.append([ligne.split('\t')[0],ligne.split('\t')[1],ligne.split('\t')[2],ligne.split('\t')[3].split('_')[0],ligne.split('\t')[5][0]])
FichUTR3.close()

# faire une boucle pour modifier le nom du gène dans lstUTR3 avec lstXRef
for i in range(0,len(lstUTR3)):
    lstUTR3[i][3]=lstXRef[lstRef.index(lstUTR3[i][3])][1]
    if (i % 1000) == 0:
        print(i*100/len(lstUTR3))

# changer l'ordre pour matcher celui des RiboSeq
Liste_chrom=["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY"]
utr3_final=[]
for c in Liste_chrom:
    for u in lstUTR3:
        if u[0]==c:
            utr3_final.append(u)

# enregistrer cette liste dans un fichier pour gagner du temps
f=open("UTR3_KG","w+")
for i in range(0,len(utr3_final)):
    f.write('\t'.join(utr3_final[i])+"\n")
f.close()

# Table UTR5 : (chromosome, debut, fin, nom usuel, sens)
lstUTR5 = []
FichUTR5 = open(NomFichUTR5)
for ligne in FichUTR5:
    lstUTR5.append([ligne.split('\t')[0],ligne.split('\t')[1],ligne.split('\t')[2],ligne.split('\t')[3].split('_')[0],ligne.split('\t')[5][0]])
FichUTR5.close()

# faire une boucle pour modifier le nom du gène dans lstUTR5 avec lstXRef
for i in range(0,len(lstUTR5)):
    lstUTR5[i][3]=lstXRef[lstRef.index(lstUTR5[i][3])][1]
    if (i % 1000) == 0:
        print(i*100/len(lstUTR5))

# changer l'ordre pour matcher celui des RiboSeq
Liste_chrom=["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY"]
utr5_final=[]
for c in Liste_chrom:
    for u in lstUTR5:
        if u[0]==c:
            utr5_final.append(u)

# enregistrer cette liste dans un fichier pour gagner du temps
f=open("UTR5_KG","w+")
for i in range(0,len(utr5_final)):
    f.write('\t'.join(utr5_final[i])+"\n")
f.close()

# création d'un fichier RB_mieux contenant : (chromosome, position, valeur)
f=open(NomFichRiboSeq)
Première_ligne = True
f2=open("RB_mieux","w+")
Chromosome = ""
for ligne in f:
    if Première_ligne:
        Première_ligne = False
    else:
        if ligne.startswith("variable"):
            Chromosome=(ligne.split('=')[1]).split(' ')[0]
            print(Chromosome)
        else:
            f2.write('\t'.join([Chromosome, ligne.split('\t')[0],ligne.split('\t')[1]]))
f.close()
f2.close()
'''

# chargement des fichier RB_mieux, UTR5_KG et UTR3_KG, et création d'un fichier RB_senescent
# operation ajoutée à la fin : moyenne

print("RB FILE : ", end=" ")
frb=open("RB_mieux")
rb=[]
for line in frb:
    rb.append([line.split('\t')[0],line.split('\t')[1],line.split('\t')[2].rstrip()])
frb.close()
print("LOADED")

print("UTR5 FILE : ", end=" ")
fu5=open("UTR5_KG")
u5=[]
for line in fu5:
    u5.append([line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[3],line.split('\t')[4].rstrip()])
fu5.close()
print("LOADED")

print("UTR3 FILE : ", end=" ")
fu3=open("UTR3_KG")
u3=[]
for line in fu3:
    u3.append([line.split('\t')[0],line.split('\t')[1],line.split('\t')[2],line.split('\t')[3],line.split('\t')[4].rstrip()])
fu3.close()
print("LOADED")

detected=False
Liste_chrom=["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY"]
RB_chrom=[] # pour donner le nombre de UTR par chromosome
for i in Liste_chrom: 
    RB_chrom.append([i,0])

'''
for c in Liste_chrom:
    print(c,": ", end=" ")
    nbmatches = 0
    # positionnement sur le début du chromosome dans les 2 tables
    iutr=0
    irb=0
    while (Liste_chrom.index(rb[irb][0])<Liste_chrom.index(c)):
        irb+=1
    while (Liste_chrom.index(u5[iutr][0])<Liste_chrom.index(c)):
        iutr+=1

    while ( iutr<len(u5) and irb< len(rb) and u5[iutr][0] == rb[irb][0] == c):
        # incrément de rb jusqu'au prochain UTR ou à la fin du chromosome
        while (irb<(len(rb)-1) and rb[irb][0]==c and int(rb[irb][1]) < int(u5[iutr][1])):
            irb+=1
        if (int(u5[iutr][1]) <= int(rb[irb][1]) <= int(u5[iutr][2])):
            nbmatches+=1
            vals=[]
            detected=True
        while (int(u5[iutr][1]) <= int(rb[irb][1]) <= int(u5[iutr][2])):
            vals.append(float(rb[irb][2]))
            irb+=1
        # calcul de la moyenne et écriture dans u5
        if detected:
            u5[iutr].append(str(mean(vals)))
            detected=False
        iutr+=1
    print(nbmatches)

# écriture de u5 dans un fichier UTR5_RB
f=open("UTR5_RB_Sen","w+")
for i in u5:
    f.write('\t'.join(i)+'\n')
f.close()
'''

# Pour UTR 3'
detected=False
for c in Liste_chrom:
    print(c,": ", end=" ")
    nbmatches = 0
    # positionnement sur le début du chromosome dans les 2 tables
    iutr=0
    irb=0
    while (Liste_chrom.index(rb[irb][0])<Liste_chrom.index(c)):
        irb+=1
    while (Liste_chrom.index(u3[iutr][0])<Liste_chrom.index(c)):
        iutr+=1

    while ( iutr<len(u3) and irb< len(rb) and u3[iutr][0] == rb[irb][0] == c):
        # incrément de rb jusqu'au prochain UTR ou à la fin du chromosome
        while (irb<(len(rb)-1) and rb[irb][0]==c and int(rb[irb][1]) < int(u3[iutr][1])):
            irb+=1
        if (int(u3[iutr][1]) <= int(rb[irb][1]) <= int(u3[iutr][2])):
            nbmatches+=1
            vals=[]
            detected=True
        while (int(u3[iutr][1]) <= int(rb[irb][1]) <= int(u3[iutr][2])):
            vals.append(float(rb[irb][2]))
            irb+=1
        # calcul de la moyenne et écriture dans u3
        if detected:
            u3[iutr].append(str(mean(vals)))
            detected=False
        iutr+=1
    print(nbmatches)

# écriture de u3 dans un fichier UTR5_RB
f=open("UTR3_RB_Sen","w+")
for i in u3:
    f.write('\t'.join(i)+'\n')
f.close()
