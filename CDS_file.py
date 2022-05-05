# Création d'un fichier contenant les chr, debut, fin, nom du gène, pour les CDS de chaque gène.

'''
Principe : création d'un tableau contenant les extrémités des CDS, le chromosome,
le brin, et finalement la moyenne et l'écart-type des valeurs Ribo-Seq
'''


# ouverture des fichiers
NomFichRiboSeqNorm = 'GSM1047589_C_Tx_alignment_reads_coverage_per_position_bothStrands_SF_0.93_.bed'
NomFichRiboSeqsen = 'GSM1047590_C.RAS_Tx_alignment_reads_coverage_per_position_bothStrands_SF_1.62_.bed'
NomFichCDS = 'CDS_hg19'
NomFichXRef = 'kgXref.txt'

# Table XRef : (nom UCSC;nom usuel)
lstXRef = []
lstRef = []
FichXRef = open(NomFichXRef)
for ligne in FichXRef:
    lstXRef.append([ligne.split('\t')[0],ligne.split('\t')[4]])
    lstRef.append(ligne.split('\t')[0])
FichXRef.close()

# Table CDS : (chromosome, debut, fin, nom UCSC, sens)
lstCDS = []
FichCDS = open(NomFichCDS)
for ligne in FichCDS:
    lstCDS.append([ligne.split('\t')[0],ligne.split('\t')[1],ligne.split('\t')[2],ligne.split('\t')[3].split('_')[0],ligne.split('\t')[5][0]])
FichCDS.close()

# faire une boucle pour modifier le nom du gène UCSC -> usuel dans lstCDS avec lstXRef
for i in range(0,len(lstCDS)):
    lstCDS[i][3]=lstXRef[lstRef.index(lstCDS[i][3])][1]

# changer l'ordre pour matcher celui des RiboSeq
Liste_chrom=["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY"]
CDS_final=[]
for c in Liste_chrom:
    for u in lstCDS:
        if u[0]==c:
            CDS_final.append(u)

# enregistrer cette liste dans un fichier pour gagner du temps
f=open("CDS_KG","w+")
for i in range(0,len(CDS_final)):
    f.write('\t'.join(CDS_final[i])+"\n")
f.close()