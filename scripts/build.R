##############
##
## This script does all the magic. Run with `Rscript scripts/build.R`
## 1-prepare build environment
## 2-statically build the root of the website
## 3-for each uniprot protein, fills in the template/protein.Rmd as <<uniprot_id>>.Rmd and renders it to html
##   3.1-the metadata section
##   3.2-pass data to template
##   3.3-save RMD file and compile
## 4-generate help like in distill:::write_search_json("staging", rmarkdown:::site_config("staging"))
## 5-rename the build directory _site to docs (where github expects the stuff to be)
## 6-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
##
###############
library(parallel)
CORES = 64

## 1-prepare build environment
unlink(c("staging", "docs"), recursive=TRUE)
dir.create("staging")
file.copy(c("_site.yml", "_footer.html", "about.Rmd", "datasets.Rmd", "downloads.Rmd", "index.Rmd", "images"),
          "staging", recursive=TRUE)

## 2-statically build the root of the website
rmarkdown::render_site("staging")#, quiet=TRUE)

## 3-for each protein listed in data/uniprot.tab
proteins <- read.delim("data/uniprot.tab")
phospho  <- read.delim(gzfile("data/Phospho_database.processed.txt.gz"))
template <- paste(readLines("templates/protein.Rmd"), collapse="\n")

i <- proteins$Gene.names...primary.. == ""
proteins[i, "Gene.names...primary.."] <- proteins$Entry.name[i]
phospho$Mod..position <- sapply(strsplit(phospho$Mod..position, ";"), function(x) paste(sort(unique(x)), collapse=";"))

## 3-fills in the template/protein.Rmd as <<uniprot_id>>.Rmd and renders it to html
system.time({
  mclapply(split(proteins, cut(seq_len(nrow(proteins)), CORES)), function(proteins) {
    # hack for parallel pandoc processing
    clean_tmpfiles_mod <- function() { invisible(0) }
    assignInNamespace("clean_tmpfiles", clean_tmpfiles_mod, ns = "rmarkdown")
    
    # copy build environment to a temp location
    dir.create(cwd <- tempfile("buildenv"))
    file.copy("staging", cwd, recursive=TRUE)
    
    # process protein RMD one by one
    invisible({
      treatments <- length(unique(phospho$Treatment))
      lapply(seq_len(nrow(proteins)), function(i) {
          # pass data to template
          j <- grepl(proteins$Entry[i], phospho$Uniprot.IDs)
          saveRDS(phospho[j, ], tmp)

          # the metadata section
          x <- gsub("<<GENE.ID>>", proteins$Gene.names...primary..[i],
               gsub("<<PROTEIN.ID>>", proteins$Entry[i],
               gsub("<<GENE.NAME>>", proteins$Gene.names[i],
               gsub("<<PROTEIN.NAME>>", proteins$Protein.names[i],
               gsub("<<FUNCTION>>", proteins$Function..CC.[i],
               gsub("<<KEYWORDS>>", proteins$Keywords[i],
               gsub("<<DATA>>", tmp <- tempfile(proteins$Entry[i]),
               gsub("<<FIG_HEATMAP_HEIGHT>>", 3 + round(1/3 * sum(j) / treatments), template))))))))

          # save RMD file and compile
          rmd  <- file.path(cwd, "staging", paste0(proteins$Entry[i], ".Rmd"))
          html <- paste0(proteins$Entry[i], ".html")
          writeLines(x, rmd)
          try(suppressWarnings(rmarkdown::render_site(rmd, envir=new.env(), quiet=TRUE)))
          file.remove(rmd, tmp)
      })
    })
    
    # merge staging folders from each core
    invisible({
      lapply(list.files(file.path(cwd, "staging/_site"), pattern="\\.html$", full=TRUE), file.copy, to="staging/_site")
      file.copy(file.path(cwd, "staging/_site/site_libs"), to="staging/_site", recursive=TRUE)
      unlink(cwd)
    })
  }, mc.cores=CORES)
})

## 4-generate help like in distill:::write_search_json("staging", rmarkdown:::site_config("staging"))
h <-  list(
  articles=data.frame(
    path       =paste0(proteins$Entry, ".html"),
    title      =proteins$Gene.names...primary..,
    description=paste(proteins$Gene.names, "❯", proteins$Protein.names)),
  collections=data.frame()
)
writeLines(jsonlite::toJSON(h), "staging/_site/search.json")
  
## 5-rename the build directory _site to docs (where github expects the stuff to be)
file.rename(file.path("staging", "_site"), "docs")
unlink("staging", recursive=TRUE)

## 6-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
