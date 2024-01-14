library(dada2)

path_to_fastqs <- "/Users/Anna/OneDrive/Cours/Master/M1/Semestre 8/stage/DADA2/dezippes"
# définir le chemin pour aller aux Fastas 

list.files(path_to_fastqs)

fnFs <- sort(list.files(path_to_fastqs,
                        pattern = "_1.fastq",
                        full.names = TRUE)) 
# définir nvlle variable qui contient tous les fichiers des fasta

# seulement sélectionner les noms de fichiers qui se finissent par ce qu'on met en pattern


fnRs <- sort(list.files(path_to_fastqs,
                        pattern = "_2.fastq",
                        full.names = TRUE))
# same pour les séquences inverses

sample.names <- basename(fnFs) |>
  strsplit(split = "_") |> # diviser la chaîne de carac selon le modèle mis entre ""
  sapply(head, 1) # appliquer une fonction à chaque éléments d'une liste
sample.names

basename(fnFs) |>
  strsplit(split = "_") |> # séparer les noms des fichiers en deux vecteurs
  head()

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

filtFs <- file.path(path_to_fastqs, "filtered", paste0(sample.names, "_1_filt.fastq"))
filtRs <- file.path(path_to_fastqs, "filtered", paste0(sample.names, "_2_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20,20), truncLen=c(240,200),
                     maxN=0, maxEE = c(2,2),truncQ=2, rm.phix=TRUE,
                     compress=FALSE, multithread=FALSE)
head(out)


#pour enlever le bruit de fond de notre data, on doit utiliser une erreur modèle, on peut l'apprendre directement par la fonction qui suit 
errF <- learnErrors(filtFs, multithread=FALSE)

errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# et là on run la fonction pour enlever le bruit de fond
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)


dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# à ce stade on a les ASV et on sait le nombre de reads de chaque échantillon, donc on peut faire une ASV table (table de comptage)

table(nchar(getSequences(seqtab)))

#les chimères sont les séquences artefacts formées par 2 ou + séquences biologiques qui sont mal jointes
# on trouve et supprime les bimères (deux chimères parents) :

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Anna/OneDrive/Cours/Master/M1/Semestre 8/stage/DADA2/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)


taxa <- addSpecies(taxa, "C:/Users/Anna/OneDrive/Cours/Master/M1/Semestre 8/stage/DADA2/silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#exporter en tant qu'objet R
export_folder <- ("export")

if (!dir.exists(export_folder)) dir.create(export_folder, recursive = TRUE)

saveRDS(object = seqtab.nochim,
        file = file.path(export_folder, "seqtab.nochim.rds"))

saveRDS(object = taxa.print,
        file = file.path(export_folder, "taxa.rds"))

write.csv(taxa.print, "tableASVfinale.csv")
