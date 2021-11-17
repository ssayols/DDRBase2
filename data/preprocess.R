library(parallel)
CORES = 4
options(mc.cores=CORES)

x <- read.delim(gzfile("data/Phospho_database.txt.gz"))
x$id <- make.names(x$Gene.names, unique=TRUE)

# filter out non-canonical isoforms
uniprot <- read.delim(gzfile("data/uniprot.tab.gz"))  # canonical == uniprot$Entry
proteins2 <- do.call(rbind, mcMap(strsplit(x$Proteins, ";"), seq_len(nrow(x)), f=function(entry, i) {
  j <- which(entry %in% uniprot$Entry)
  c(Proteins =paste(entry[j], collapse=";"),
    Positions=paste(unlist(strsplit(x$Positions.within.proteins[i], ";"))[j], collapse=";"))
}))
x$Proteins <- proteins2[, "Proteins"]
x$Positions.within.proteins <- proteins2[, "Positions"]

# convert wide to long
x <- reshape(x,
             direction="long", 
             varying  =list(grep("_logFC$", colnames(x), value=TRUE), grep("_adj.P.Val$", colnames(x), value=TRUE)),
             timevar  ="Treatment",
             times    =c("MMS", "HU", "H2O2", "FA", "ETO", "UV", "APH", "CPT", "Xray", "Gem", "ArO2"),
             v.names  =c("logFC", "adj.P.Val"),
             idvar    ="id")

# Additional information on the experiment/cell_line/treatment
x$Enrichment <- "phospho"
x$Cell_line  <- "U2OS"

specifics <- c("APH" ="10 µM, 2h",
               "ArO2"="140 µM, 2h",
               "CPT" ="21 µM, 2h",
               "ETO" ="26 µM, 2h",
               "FA"  ="200 µM, 2h",
               "Gem" ="10 µM, 2h",
               "H2O2"="400 µM, 2h",
               "HU"  ="2 mM, 2h",
               "MMS" ="0.0044%, 2h",
               "UV"  ="12 J/m2, 1h recovery",
               "Xray"="2.4 Gy, 2h recovery")
x$Specifics <- specifics[match(x$Treatment, names(specifics))]

treatments <- c("APH" ="Aphidicolin",
                "ArO2"="Sodium Arsenite",
                "CPT" ="Camptothecin",
                "ETO" ="Etoposide",
                "FA"  ="Formaldehyd",
                "Gem" ="Gemcitabine",
                "H2O2"="Hydrogen peroxide",
                "HU"  ="Hydroxyurea",
                "MMS" ="Methyl methanesulfonate",
                "UV"  ="UV-C",
                "Xray"="X-ray")
x$Treatment <- factor(x$Treatment, levels=names(treatments), labels=treatments)

# rename output columns to something nicer
dcols <- c(`Gene name`         ="Gene.names",
           `Uniprot IDs`       ="Proteins",
           `Protein name`      ="Protein.names",
           `Enrichment`        ="Enrichment",
           `Treatment`         ="Treatment",
           `Specifics`         ="Specifics",
           `Cell line`         ="Cell_line",
           `Mod. position`     ="Positions.within.proteins",
           `Amino Acid`        ="Amino.acid",
           `Sequence window`   ="Sequence.window",
           `Localization prob.`="Phospho..STY..Probabilities",
           `Multiplicity`      ="multiplicity",
           `log2FC`            ="logFC",
           `FDR`               ="adj.P.Val")
x <- x[, dcols]
colnames(x) <- names(dcols)

# save processed table
write.table(x, gzfile("data/Phospho_database.processed.txt.gz"), row.names=FALSE, sep="\t")
