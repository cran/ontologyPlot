## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev="svg", fig.width=7, fig.height=5, dev="svg", fig.align="center")

## ----fig.width=6.15, fig.height=5.67-------------------------------------
library(ontologyIndex)
library(ontologyPlot)
data(go)

genes <- list(
	A0A087WUJ7=c("GO:0004553", "GO:0005975"),
	CTAGE8=c("GO:0016021"),
	IFRD2=c("GO:0003674", "GO:0005515", "GO:0005634"),
	OTOR=c("GO:0001502", "GO:0005576", "GO:0007605"),
	TAMM41=c("GO:0004605", "GO:0016024", "GO:0031314", "GO:0032049"),
	ZZEF1=c("GO:0005509", "GO:0008270")
)

plot_annotation_grid(go, term_sets=genes)

