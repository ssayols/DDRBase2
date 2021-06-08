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

## prepare build environment
unlink(c("staging", "docs"), recursive=TRUE)
dir.create("staging")
dir.create("staging/proteins")  # where the rendered proteins go
file.copy(c("_site.yml", "_footer.html", "about.Rmd", "datasets.Rmd", "downloads.Rmd", "index.Rmd", "images"),
          "staging", recursive=TRUE)

##  1-statically build the root of the website
rmarkdown::render_site("staging")#, quiet=TRUE)

##  1-for each protein listed in data/uniprot.tab
proteins <- head(read.delim("data/uniprot.tab"))
template <- paste(readLines("templates/protein.Rmd"), collapse="\n")

##   1.1-fills in the template/protein.Rmd as <<uniprot_id>>.Rmd and renders it to html
invisible({
  mclapply(seq_len(nrow(proteins)), function(i) {
    # the metadata section
    x <- gsub("<<GENE.ID>>", proteins$Gene.names...primary..[i],
         gsub("<<PROTEIN.ID>>", proteins$Entry[i],
         gsub("<<GENE.NAME>>", proteins$Gene.names[i],
         gsub("<<PROTEIN.NAME>>", proteins$Protein.names[i],
         gsub("<<FUNCTION>>", proteins$Function..CC.[i],
         gsub("<<KEYWORDS>>", proteins$Keywords[i], template))))))
    
    # test plot (as base64)
    png(tf1 <- tempfile(fileext=".png")); plot(0); dev.off()
    txt     <- RCurl::base64Encode(readBin(tf1, "raw", file.info(tf1)[1, "size"]), "txt")
    html    <- sprintf('<img src="data:image/png;base64,%s">', txt)
    x       <- gsub("<<TEST.PLOT>>", html, x)
    
    # save file
    rmd <- file.path("staging", paste0(proteins$Entry[i], ".Rmd"))
    writeLines(x, rmd)
    rmarkdown::render_site(rmd)#, quiet=TRUE)
    file.copy(file.path("staging/_site", paste0(proteins$Entry[i], ".html")), "staging/_site/proteins")
    file.remove(rmd)
  }, mc.cores=1)
})

## generate help like in distill:::write_search_json("staging", rmarkdown:::site_config("staging"))
h <-  data.frame(title      =proteins$Gene.names...primary..,
                 description=paste(proteins$Gene.names, "❯", proteins$Protein.names))
h$path <- file.path("proteins", paste0(proteins$Entry, ".html"))
writeLines(jsonlite::toJSON(h), "staging/_site/search.json")
  
## 
file.rename(file.path("staging", "_site"), "docs")
unlink("staging", recursive=TRUE)

##  3-minify the search, to do it only on the description of the page
#system("sed -i '/^      \"contents\": /d;/^      \"author\": /d;/^      \"last_modified\": /d' docs/search.json")

##  3-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
