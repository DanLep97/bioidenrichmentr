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
enrichedGOdata = data.frame(
  bait = character(0),
  goID = character(0),
  term = character(0),
  score = character(0),
  ont = character(0),
  stringsAsFactors = FALSE
)
enrichedPCdata = data.frame(
  bait = character(0),
  pcID = character(0),
  pcCount = integer(0),
  complexName = character(0),
  uniprotID = character(0),
  stringsAsFactors = FALSE
)

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
graphs = vector("list", length=length(baits))
names(graphs) = baits

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
filterGenes = function(GOids, ont, bait, goData) {
  
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
    
    
    ## FIND UNIQUE GENES IN TERM
    genesInGoTerm = gene2GOdf[which(
      gene2GOdf$goID == goTerm &
        gene2GOdf$uniprotID %in% annotatedGenesInGoTerm
    ), "uniprotID"]
    
    ## FIND REDENDANCIES IN ANCESTORS AND CHILDREN
    ancestors = unname(totAncestors[[goTerm]])
    offspring = unname(totOffspring[[goTerm]])
    
    offspringGenes = gene2GOdf[which(
      gene2GOdf$goID %in% offspring &
        gene2GOdf$uniprotID %in% annotatedGenesInGoTerm
    ), "uniprotID"]
    offspringGenes = unique(offspringGenes)
    
    ancestorsGenes = gene2GOdf[which(
      gene2GOdf$goID %in% ancestors &
        gene2GOdf$uniprotID %in% annotatedGenesInGoTerm
    ), "uniprotID"]
    ancestorsGenes = unique(ancestorsGenes)
    
    for(j in 1:length(annotatedGenesInGoTerm)) {
      gene = annotatedGenesInGoTerm[j]
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
      genesOfTerm <<- rbind(genesOfTerm, data.frame(
        goID = goTerm,
        gene = gene,
        uniqueness = str
      ))
    }
  }
}
enrichGOterms = function(bait, ont) {
  minNumberOfGenesPerTerms = 5
  myInterestingGenes = geneNames[[bait]]
  print(c("number of genes for selected bait", length(myInterestingGenes)))
  dataGeneNames = unique(sData$UniprotID)
  geneList = factor(as.integer(dataGeneNames %in% myInterestingGenes))
  names(geneList) = dataGeneNames
  
  GOdataBP = new("topGOdata",
                 description="topGO object for a given bait",
                 ontology=ont,
                 allGenes = geneList,
                 nodeSize= minNumberOfGenesPerTerms,
                 annot= annFUN.gene2GO,
                 gene2GO = gene2GO
  )
  weight01FisherResult = runTest(GOdataBP, statistic = "fisher") # default algorithm = weight01
  
  if (ont=="BP") {
    #add graph to the list of graphs
    dag = showSigOfNodes(GOdataBP, score(weight01FisherResult), firstSigNodes = 5, useInfo = "def")
    dotFile = tempfile()
    toFile(dag$complete.dag, filename = dotFile)
    dotStr = readLines(dotFile)
    graphs[[bait]] <<- dotStr
    unlink(dotFile)
  }
  
  allRes = GenTable(
    GOdataBP,
    weight = weight01FisherResult,
    topNodes = 20,
    orderBy="weight"
  )
  filterGenes(allRes$GO.ID, ont, bait, GOdataBP)
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
      stringsAsFactors = FALSE
    )
    return(returnedTable)
  }
}
enrichBaits = function(length) {
  for(i in 1:length) {
    BPterms = enrichGOterms(baits[i], "BP") # enrich biological processes
    CCterms = enrichGOterms(baits[i], "CC") # enrich cellular locations
    PCterms = enrichPC(baits[i]) # enrich protein complexes
    
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
      stringsAsFactors = FALSE)
    
    enrichedGOdata <<- rbind(enrichedGOdata, CCdf, BPdf)
    enrichedPCdata <<- rbind(enrichedPCdata, PCdf)
  }
}


#API FUNCTIONS:
#* @post fingerprint
#* @param length
function(length) {
  enrichBaits(length)
  superEnrichedGOdata = data.frame()
  superEnrichedPCdata = data.frame()
  for(i in 1:length) {
    bait = baits[i]
    
    baitGOensembl = enrichedGOdata[which(enrichedGOdata$bait == bait),]
    substrGOensembl = enrichedGOdata[which(enrichedGOdata$bait != bait),]
    
    baitPCensembl = enrichedPCdata[which(enrichedPCdata$bait == bait),]
    substrPCensembl = enrichedPCdata[which(enrichedPCdata$bait != bait),]
    
    intersectGO = na.omit(intersect(baitGOensembl$goID, substrGOensembl$goID))
    intersectPC = na.omit(intersect(baitPCensembl$pcID, substrPCensembl$pcID))
    
    for(j in 1:length(intersectGO)) {
      ind = which(
        baitGOensembl$goID == intersectGO[j]
      ) 
      baitGOensembl = baitGOensembl[-ind,]
    }
    
    for(j in 1:length(intersectPC)) {
      ind = which(
        baitPCensembl$pcID == intersectPC[j]
      ) 
      baitPCensembl = baitPCensembl[-ind,]
    }
    
  superEnrichedGOdata = rbind(superEnrichedGOdata, baitGOensembl)
  superEnrichedPCdata = rbind(superEnrichedPCdata, baitPCensembl)
  }
  
  #eliminate duplicates genes NOT ACTIVE
  indexsToRem = which(
    genesOfTerm$uniqueness == "|offspring|ancestor|inTerm" |
      duplicated(genesOfTerm$gene)
  )
  genesOfTerm[-indexsToRem,]
  
  #output the result as a list
  res = list(
    #"enrichedGenes" = enrichedGenes,
    "enrichedGenes" = genesOfTerm,
    "enrichedGOdata" = superEnrichedGOdata,
    "enrichedPCdata" = superEnrichedPCdata,
    "baits" = baits[1:length],
    "graphs" = graphs
  )
  write_json(res, "~/Desktop/sitestage/enrichedData.json")
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
