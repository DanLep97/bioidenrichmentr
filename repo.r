library(topGO)
library(readxl)
library(dplyr)
library(hablar)
library(naniar)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(gridExtra)
library(Rgraphviz)
library(threadr)

# data prep:
sData = read_excel("sData.xlsx")
gene2GO = readMappings(file="./gene2GO.map")
gene2GOdf = read.csv2("./gene2GOdf.csv")

#Enriched genes associated to GO term and its bait
genesOfTerm = data.frame()

#enrich data for the fingerprint
enrichedGOdata = data.frame()
enrichedPCdata = data.frame()

sData = sData %>%
  convert(dbl(np,
              nq,
              nfc,
              cp,
              cq,
              cfc,
              polyICp,
              polyICq,
              polyICFC,
              polyICFCB,
  ))
sData %>%
  replace_with_na(replace = list(x = c("na", "ns")))
sData$UniprotID = str_replace_all(sData$UniprotID, "-[0-9]*", "")

# each bait is populated with its preys:
baits = unique(sData$Bait)
geneNames = vector("list", length=length(baits))
names(geneNames) = baits
for (i in 1:length(baits)) {
  rows = sData %>% filter(Bait == baits[i])
  geneNames[[i]] = rows$UniprotID
}

#Graph containing the dotlist strings
graphsBP = vector("list", length=length(baits))
graphsCC = vector("list", length=length(baits))
names(graphsCC) = baits
names(graphsBP) = baits

#CORUM data prep:
corumData = read.delim2("./CORUM.txt", sep="\t", stringsAsFactors = FALSE)

#API config:
#' @filter cors
cors <- function(req, res) {
  
  res$setHeader("Access-Control-Allow-Origin", "*")
  
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200 
    return(list())
  } else {
    plumber::forward()
  }
  
}

#REPO FUNCTIONS
filterGenes = function(GOids, ont, bait, goData, geneRange) {
  
  if (ont == "CC") {
    totAncestors = as.list(GOCCPARENTS)
    totOffspring = as.list(GOCCOFFSPRING)
  } else {
    totAncestors = as.list(GOBPPARENTS)
    totOffspring = as.list(GOBPOFFSPRING)
  }
  for(i in 1:length(GOids)) {
    goTerm = GOids[i]
    annotatedGenesInGoTerm = genesInTerm(goData, goTerm)[[1]]
    sigGenesInTerm = annotatedGenesInGoTerm[annotatedGenesInGoTerm %in% geneRange]
    
    ## FIND UNIQUE GENES IN TERM
    genesInGoTerm = gene2GOdf[which(
      gene2GOdf$goID == goTerm &
        gene2GOdf$uniprotID %in% sigGenesInTerm
    ), "uniprotID"]
    
    ## FIND REDENDANCIES IN ANCESTORS AND CHILDREN
    ancestors = unname(totAncestors[[goTerm]])
    offspring = unname(totOffspring[[goTerm]])
    
    offspringGenes = gene2GOdf[which(
      gene2GOdf$goID %in% offspring &
        gene2GOdf$uniprotID %in% sigGenesInTerm
    ), "uniprotID"]
    
    ancestorsGenes = gene2GOdf[which(
      gene2GOdf$goID %in% ancestors &
        gene2GOdf$uniprotID %in% sigGenesInTerm
    ), "uniprotID"]
    ancestorsGenes = unique(ancestorsGenes)
    
    for(j in 1:length(sigGenesInTerm)) {
      gene = sigGenesInTerm[j]
      str = ""
      if (gene %in% offspringGenes) {
        str = paste(str, "offspring", sep = "|")
      }
      if (gene %in% ancestorsGenes) {
        str = paste(str, "ancestor", sep = "|")
      }
      if (gene %in% genesInGoTerm) {
        str = paste(str, "inTerm", sep = "|")
      }
      geneSymbol = sData[which(sData$UniprotID == gene),"GeneName"]
      genesOfTerm <<- rbind(genesOfTerm, data.frame(
        stringsAsFactors = FALSE,
        goID = goTerm,
        uniprotID = gene,
        geneSymbol = geneSymbol[1,"GeneName"],
        uniqueness = str
      ))
    }
  }
}
enrichGOterms = function(bait, ont,nodesOnGraph) {
  myInterestingGenes = geneNames[[bait]]
  print(c("number of genes for selected bait", length(myInterestingGenes)))
  dataGeneNames = unique(sData$UniprotID)
  geneList = factor(as.integer(dataGeneNames %in% myInterestingGenes))
  names(geneList) = dataGeneNames
  
  GOdataBP = new("topGOdata",
    description="topGO object for a given bait",
    ontology=ont,
    allGenes = geneList,
    nodeSize= 5,
    annot= annFUN.gene2GO,
    gene2GO = gene2GO
  )
  weight01FisherResult = runTest(GOdataBP, statistic = "fisher") # default algorithm = weight01
  
  #add graph to the list of graphs
  dag = showSigOfNodes(
    GOdataBP, 
    score(weight01FisherResult), 
    useFullNames = TRUE,
    useInfo = "all",
    .NO.CHAR = 60,
    firstSigNodes = nodesOnGraph
  )
  dotFile = tempfile()
  toFile(dag$complete.dag, filename = dotFile)
  dotStr = readLines(dotFile)
  if (ont == "BP") {
    graphsBP[[bait]] <<- dotStr
  } else {
    graphsCC[[bait]] <<- dotStr
  }
  unlink(dotFile)
  
  allRes = GenTable(
    GOdataBP,
    weight = weight01FisherResult,
    topNodes = 20,
    orderBy="weight"
  )
  filterGenes(allRes$GO.ID, ont, bait, GOdataBP, myInterestingGenes)
  return(allRes)
}
enrichPC = function(bait) {
  preys = geneNames[[bait]]
  preyLength = length(preys)
  CORUMdf = data.frame(
    Prey=character(0), Complex=character(0), stringsAsFactors = FALSE)
  
  for(i in 1:preyLength) {
    prey = preys[i]
    #print(prey)
    complexes = corumData[which(str_detect(corumData$subunits.UniProt.IDs., prey)),"ComplexID"]
    #print(length(complexes))
    if (length(complexes) == 0) {
      CORUMdf[nrow(CORUMdf)+1, c("Prey", "Complex")] = c(prey,NA)
    } else {
      for(j in 1:length(complexes)) {
        CORUMdf[nrow(CORUMdf)+1, c("Prey", "Complex")] = c(prey,complexes[j])
      }
    }
  }
  
  # make the occurence count table with added columns of interest as the presentation table. 
  idCounts = as.data.frame(table(CORUMdf$Complex), stringsAsFactors = FALSE)
  if(nrow(idCounts) > 0) {
    names(idCounts) = c("id", "occurences")
    #returnedTable = merge(x=corumData, y=idCounts, by.y="id", by.x="ComplexID", all.y=TRUE)
    CORUMdf = merge(x=CORUMdf, y=idCounts, by.y="id", by.x="Complex", all.y=TRUE)
    returnedTable = merge(x=corumData, y=CORUMdf, by.y="Complex", by.x="ComplexID", all.y=TRUE)
  } else {
    returnedTable = data.frame(
      ComplexName = NA,
      occurences = NA,
      ComplexID = NA,
      GeneName = NA,
      stringsAsFactors = FALSE
    )
  }
  return(returnedTable)
}
enrichBaits = function(length,nodeLength) {
  for(i in 1:length) {
    print(paste("enriching bait", baits[i], ".", i, "out of", length, sep=" "))
    
    BPterms = enrichGOterms(baits[i], "BP",nodeLength) # enrich biological processes
    CCterms = enrichGOterms(baits[i], "CC",nodeLength) # enrich cellular locations
    PCterms = enrichPC(baits[i]) # enrich protein complexes
    
    PCgeneName = data.frame()
    for (j in 1:nrow(PCterms)) {
      val = unique(sData[which(sData$UniprotID == PCterms[j, "Prey"]),"GeneName"])
      PCgeneName = rbind(PCgeneName, data.frame(val))
    }
    
    CCdf = data.frame(
      bait=baits[i], 
      goID = CCterms$GO.ID, 
      term = Term(CCterms$GO.ID),
      ont = "CC",
      score = CCterms$weight,
      stringsAsFactors = FALSE)
    BPdf = data.frame(
      bait=baits[i], 
      goID=BPterms$GO.ID, 
      term = Term(BPterms$GO.ID),
      ont = "BP",
      score = BPterms$weight,
      stringsAsFactors = FALSE)
    PCdf = data.frame(
      bait=baits[i], 
      pcID=PCterms$ComplexID, 
      pcCount=PCterms$occurences,
      complexName = PCterms$ComplexName,
      uniprotID = PCterms$Prey,
      geneName = PCgeneName$GeneName,
      stringsAsFactors = FALSE)
    
    enrichedGOdata <<- rbind(enrichedGOdata, CCdf, BPdf)
    enrichedPCdata <<- rbind(enrichedPCdata, PCdf)
  }
}


#API FUNCTIONS:
#* @post fingerprint
#* @param length
#* @param nodeLength
function(length,nodeLength) {
  enrichBaits(length,nodeLength)
  superEnrichedGOdata = data.frame()
  for(i in 1:length) {
    bait = baits[i]
    
    baitGOensembl = enrichedGOdata[which(enrichedGOdata$bait == bait),]
    substrGOensembl = enrichedGOdata[which(enrichedGOdata$bait != bait),]
    
    intersectGO = na.omit(intersect(baitGOensembl$goID, substrGOensembl$goID))
    
    for(j in 1:length(intersectGO)) {
      ind = which(
        baitGOensembl$goID == intersectGO[j]
      ) 
      baitGOensembl = baitGOensembl[-ind,]
    }
    
  superEnrichedGOdata = rbind(superEnrichedGOdata, baitGOensembl)
  }
  
  #enrich GO and PC genes
  duplicatesInEnrichedGOgenes = which(duplicated(genesOfTerm$uniprotID) & 
                                        genesOfTerm$goID %in% superEnrichedGOdata$goID
  )
  duplicatesInEnrichedPCgenes = which(duplicated(enrichedPCdata$uniprotID))
  
  GOduplicates = genesOfTerm[duplicatesInEnrichedGOgenes, "uniprotID"]
  PCduplicates = enrichedPCdata[duplicatesInEnrichedPCgenes, "uniprotID"]
  
  for (i in 1:nrow(genesOfTerm)) { #for genes of GO
    if (genesOfTerm[i, "uniprotID"] %in% GOduplicates) {
      genesOfTerm[i, "uniqueness"] = paste(genesOfTerm[i, "uniqueness"], "duplicated", sep="|")
    }
  }
  for (i in 1:nrow(enrichedPCdata)) { #for genes of PC
    if (enrichedPCdata[i, "uniprotID"] %in% PCduplicates) {
      enrichedPCdata[i, "uniqueness"] = "duplicated"
    } else {
      enrichedPCdata[i, "uniqueness"] = "unique"
    }
  }
  
  #output the result as a list
  res = list(
    #"enrichedGenes" = enrichedGenes,
    "enrichedGenes" = genesOfTerm,
    "enrichedGOdata" = superEnrichedGOdata,
    "enrichedPCdata" = enrichedPCdata,
    "baits" = baits[1:length],
    "graphsBP" = graphsBP,
    "graphsCC" = graphsCC
  )
  write_json(res, "~/Desktop/sitestage/enrichedData.json")
  genesOfTerm = data.frame()
  superEnrichedGOdata = data.frame()
  enrichedPCdata = data.frame()
  graphsBP = vector("list", length=length(baits))
  graphsCC = vector("list", length=length(baits))
  names(graphsCC) = baits
  names(graphsBP) = baits
  return(res)
}

#* @get corum
#* @param bait used for protein complex ordering
function(bait) {
  preys = geneNames[[bait]]
  preyLength = length(preys)
  CORUMdf = data.frame(Prey=character(0), Complexes=character(0), stringsAsFactors = FALSE)

  for(i in 1:preyLength) {
    #prey = CORUMdf[i,"Prey"]
    prey = str_replace_all(preys[i], "-[0-9]*", "")
    #print(prey)
    complexes = corumData[which(str_detect(corumData$subunits.UniProt.IDs., prey)),"ComplexID"]
    #print(length(complexes))
    if (length(complexes) == 0) {
      CORUMdf[nrow(CORUMdf)+1, c("Prey", "Complexes")] = c(prey,NA)
    } else {
      for(j in 1:length(complexes)) {
        CORUMdf[nrow(CORUMdf)+1, c("Prey", "Complexes")] = c(prey,complexes[j])
      }
    }
  }
  
  # make the occurence count table with added columns of interest as the presentation table. 
  idCounts = as.data.frame(table(CORUMdf$Complexes), stringsAsFactors = FALSE)
  names(idCounts) = c("id", "occurences")
  idCounts = idCounts[which(idCounts$occurences>=2),]
  returnedTable = merge(x=corumData, y=idCounts, by.y="id", by.x="ComplexID", all.y=TRUE)
  return(returnedTable)
}

#* @get gopb
#* @param bait used for BP and CC
function(bait) {
  myInterestingGenes = geneNames[[bait]]
  print(c("number of genes for selected bait", length(myInterestingGenes)))
  dataGeneNames = unique(sData$UniprotID)
  geneList = factor(as.integer(dataGeneNames %in% myInterestingGenes))
  names(geneList) = dataGeneNames
  
  GOdataBP = new("topGOdata",
                 description="topGO object for a given bait",
                 ontology="BP",
                 allGenes = geneList,
                 nodeSize= minNumberOfGenesPerTerms,
                 annot= annFUN.gene2GO,
                 gene2GO = gene2GO
  )
  weight01FisherResult = runTest(GOdataBP, statistic = "fisher") # default algorithm = weight01
  
  dotFile = tempfile()
  toDot(graph(GOdataBP), dotFile)
  dotStr = readLines(dotFile);
  return(dotStr)
}
