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
unlink(c("_site", "docs"), recursive=TRUE)

##  1-for each protein listed in data/uniprot.tab
proteins <- head(read.delim("data/uniprot.tab"))
template <- paste(readLines("templates/protein.Rmd"), collapse="\n")

##   1.1-fills in the template/protein.Rmd as <<uniprot_id>>.Rmd
for(i in seq_len(nrow(proteins))) {
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
  writeLines(x, paste0(proteins$Entry[i], ".Rmd"))
}

##  1-statically build the root of the website
rmarkdown::render_site()
# TODO: remove protein Rmd files

## 
file.rename("_site", "docs")
unlink(c(file.path("docs", "data"),  file.path("docs", "scripts"),  file.path("docs", "templates")), recursive=TRUE)

##  3-minify the search, to do it only on the description of the page
system("sed -i '/^      \"contents\": /d;/^      \"author\": /d;/^      \"last_modified\": /d' docs/search.json")

##  3-`git commit -am'rebuild' && git push`, and go to http://ssayols.github.io/DDRBase2
