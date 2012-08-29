#!/usr/bin/python
import os, sys
import rpy


def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " universe-topGO.table  testset.ids"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inUniverse, inTestset = sys.argv[1:3]
  return inUniverse, inTestset


def init_R():
  R = rpy.r
  try:
    R.library('topGO')
  except:
    try: 
      R.source("http://bioconductor.org/biocLite.R")
      R.biocLite('topGO')
      R.library('topGO')
    except:
      print "Problem importing R libraries."
      sys.exit()

  R('if(!isGeneric("GOFisherTestUnder")) setGeneric("GOFisherTestUnder", function(object) standardGeneric("GOFisherTestUnder"))')
  R('setMethod("GOFisherTestUnder", "classicCount", function(object) { contMat <- contTable(object); if(all(contMat == 0)) p.value <- 1 else p.value <- fisher.test(contMat, alternative = "less")$p.value; return(p.value) })')
  return R


def main():
  inUniverse, inTestset = plausi()
  R = init_R()
  R('GOmap = readMappings(file = "' + inUniverse + '")')
  R('refset  = names(GOmap)')
  R('testset = scan(file="' + inTestset + '", what=character())')
  R('genes_of_interest = factor(as.integer(refset %in% testset))')
  R('names(genes_of_interest) <- refset')
  for ontology in ["MF", "BP", "CC"]:
    R('tgData = new("topGOdata", ontology = "' + ontology + '", allGenes = genes_of_interest, annot = annFUN.gene2GO, gene2GO = GOmap)')
    R('fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")')
    R('fisherResCor = p.adjust(score(fisherRes), method="fdr")')
    R('weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")')
    R('weightResCor = p.adjust(score(weightRes), method="fdr")')
    R('allRes    = GenTable(tgData, classic=fisherRes, weight=weightRes, orderBy="weight", ranksOf="classic", topNodes=150)')
    R('allRes$fisher.FDR = fisherResCor[allRes$GO.ID]')
    R('allRes$weight.FDR = weightResCor[allRes$GO.ID]')
    R('write.csv(allRes, "topGO.over.Sig.' + ontology + '.csv")')

    R('tgData = new("topGOdata", ontology = "' + ontology + '", allGenes = genes_of_interest, annot = annFUN.gene2GO, gene2GO = GOmap)')
    R('test.stat <- new("classicCount", testStatistic = GOFisherTestUnder, name ="Fisher test underrepresentation")')
    R('fisherRes <- getSigGroups(tgData, test.stat)')
    R('fisherResCor = p.adjust(score(fisherRes), method="fdr")')
    R('test.stat <- new("weightCount", testStatistic = GOFisherTestUnder, name ="Fisher test underrepresentation")')
    R('weightRes <- getSigGroups(tgData, test.stat)')
    R('weightResCor = p.adjust(score(weightRes), method="fdr")')
    R('allRes    = GenTable(tgData, classic=fisherRes, weight=weightRes, orderBy="weight", ranksOf="classic", topNodes=150)')
    R('allRes$fisher.FDR = fisherResCor[allRes$GO.ID]')
    R('allRes$weight.FDR = weightResCor[allRes$GO.ID]')
    R('write.csv(allRes, "topGO.under.Sig.' + ontology + '.csv")')
 
main()
