##############
##
## This script does all the magic. Run with `Rscript scripts/build.R`
##  2-for each protein listed in data/uniprot.tab
##   2.1-fills in the template/protein.Rmd as <<uniprot_id>>.Rmd
##   2.2-builds the resulting Rmd files
##   2.2-copies the resulting html files to docs/proteins
##  1-statically build the root of the website
##    `rmarkdown::render_site()`
##    `file.rename("_site", "docs")`   # since `docs` is where github.io serves the site
##  3-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
##
###############
library(parallel)
CORES = 8

## 1-prepare build environment
unlink(c("staging", "docs"), recursive=TRUE)
dir.create("staging")
file.copy(c("_site.yml", "_footer.html", "about.Rmd", "datasets.Rmd", "downloads.Rmd", "index.Rmd", "images"),
          "staging", recursive=TRUE)

## 2-statically build the root of the website
rmarkdown::render_site("staging")#, quiet=TRUE)

## 3-for each protein listed in data/uniprot.tab
proteins <- head(read.delim("data/uniprot.tab"))
phospho  <- read.delim(gzfile("data/Phospho_database.processed.txt.gz"))
template <- paste(readLines("templates/protein.Rmd"), collapse="\n")

## 3.1-fills in the template/protein.Rmd as <<uniprot_id>>.Rmd and renders it to html
dir.create("staging/_site/proteins")
system.time({
#  lapply(split(proteins, cut(seq_len(nrow(proteins)), CORES)), function(proteins) {
#    # copy build environment to a temp location
#    dir.create(cwd <- tempfile())
#    file.copy("staging", cwd, recursive=TRUE)
    
    cwd="."
    # process protein RMD one by one
    lapply(seq_len(nrow(proteins)), function(i) {
      # the metadata section
      x <- gsub("<<GENE.ID>>", proteins$Gene.names...primary..[i],
           gsub("<<PROTEIN.ID>>", proteins$Entry[i],
           gsub("<<GENE.NAME>>", proteins$Gene.names[i],
           gsub("<<PROTEIN.NAME>>", proteins$Protein.names[i],
           gsub("<<FUNCTION>>", proteins$Function..CC.[i],
           gsub("<<KEYWORDS>>", proteins$Keywords[i],
           gsub("<<DATA>>", tmp <- tempfile(proteins$Entry[i]), template)))))))
      
      # pass data to template
      saveRDS(phospho[grepl(proteins$Entry[i], phospho$Uniprot.IDs), ], tmp)
      
      # save RMD file and compile
      rmd  <- file.path(cwd, "staging", paste0(proteins$Entry[i], ".Rmd"))
      html <- paste0(proteins$Entry[i], ".html")
      writeLines(x, rmd)
      rmarkdown::render_site(rmd, quiet=TRUE)
      file.rename(file.path(cwd, "staging/_site", html), file.path(cwd, "staging/_site/proteins", html))
      file.remove(rmd)
    })
    
    # merge staging folders from each core
    #file.copy(file.path(cwd, "staging/_site/proteins"), "staging/_site/proteins")
#  }, mc.cores=CORES)
})

## 4-generate help like in distill:::write_search_json("staging", rmarkdown:::site_config("staging"))
h <-  list(
  articles=data.frame(
    path       =file.path("proteins", paste0(proteins$Entry, ".html")),
    title      =proteins$Gene.names...primary..,
    description=paste(proteins$Gene.names, "❯", proteins$Protein.names)),
  collections=data.frame()
)
writeLines(jsonlite::toJSON(h), "staging/_site/search.json")
  
## 5-rename the build directory _site to docs (where github expects the stuff to be)
file.rename(file.path("staging", "_site"), "docs")
unlink("staging", recursive=TRUE)

## 6-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
