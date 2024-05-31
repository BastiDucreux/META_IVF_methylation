library(readxl)
library(GenomicRanges)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(bumphunter)
library("TxDb.Hsapiens.UCSC.hg18.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(dplyr)
library(rtracklayer)

set.seed(100)

#DMP analysis

YeungDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Yeung DMP"))
BarberetDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Barberet DMP"))
DucreuxDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Ducreux DMP"))
TobiDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Tobi DMP"))
NovakovicneoDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Novakovic DMP neonates"))
MelamedDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Melamed DMP"))
KatariDMPcord <- data.frame(read_excel("DMP annot.xlsx",sheet = "Katari DMP Cord blood"))
KatariDMPplacenta <- data.frame(read_excel("DMP annot.xlsx",sheet = "Katari DMP placenta"))
ElHajjDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "El Hajj DMP"))
CamaraschiDMP <- data.frame(read_excel("DMP annot.xlsx",sheet = "Camaraschi DMP"))

ListDMP=list("YeungDMP"=YeungDMP,"BarberetDMP"=BarberetDMP,"DucreuxDMP"=DucreuxDMP,
             "TobiDMP"=TobiDMP,"NovakovicneoDMP"=NovakovicneoDMP,"MelamedDMP"=MelamedDMP,"KatariDMPcord"=KatariDMPcord,"KatariDMPplacenta"=KatariDMPplacenta,"ElHajjDMP"=ElHajjDMP,"CamaraschiDMP"=CamaraschiDMP)

for (k in 1:length(ListDMP)){
  for (j in k:length(ListDMP)){
    if (names(ListDMP[k])!=names(ListDMP[j])){
      cat(c(names(ListDMP[k]),names(ListDMP[j])))
      cat(intersect(ListDMP[[k]]$CpG,ListDMP[[j]]$CpG),sep = '\n')
      cat(sep = '\n')
    }
  }
} 

genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

YeungDMP$Gene=matchGenes(YeungDMP,genes)$name
BarberetDMP$Gene=matchGenes(BarberetDMP,genes)$name
DucreuxDMP$Gene=matchGenes(DucreuxDMP,genes)$name
TobiDMP$Gene=matchGenes(TobiDMP,genes)$name
NovakovicneoDMP$Gene=matchGenes(NovakovicneoDMP,genes)$name
MelamedDMP$Gene=matchGenes(MelamedDMP,genes)$name
ElHajjDMP$Gene=matchGenes(ElHajjDMP,genes)$name
CamaraschiDMP$Gene=matchGenes(CamaraschiDMP,genes)$name

genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg18.knownGene)

KatariDMPcord$Gene=matchGenes(KatariDMPcord,genes)$name
KatariDMPplacenta$Gene=matchGenes(KatariDMPplacenta,genes)$name

#selection genes with >=2 DMPs to match Katari et al. criteria

YeungDMP <- subset(YeungDMP,YeungDMP$Gene%in%names(which(table(YeungDMP$Gene)>=2)))
BarberetDMP <- subset(BarberetDMP,BarberetDMP$Gene%in%names(which(table(BarberetDMP$Gene)>=2)))
DucreuxDMP <- subset(DucreuxDMP,DucreuxDMP$Gene%in%names(which(table(DucreuxDMP$Gene)>=2)))
TobiDMP <- subset(TobiDMP,TobiDMP$Gene%in%names(which(table(TobiDMP$Gene)>=2)))
NovakovicneoDMP <- subset(NovakovicneoDMP,NovakovicneoDMP$Gene%in%names(which(table(NovakovicneoDMP$Gene)>=2)))
MelamedDMP <- subset(MelamedDMP,MelamedDMP$Gene%in%names(which(table(MelamedDMP$Gene)>=2)))
KatariDMPcord <- subset(KatariDMPcord,KatariDMPcord$Gene%in%names(which(table(KatariDMPcord$Gene)>=2)))
KatariDMPplacenta <- subset(KatariDMPplacenta,KatariDMPplacenta$Gene%in%names(which(table(KatariDMPplacenta$Gene)>=2)))
ElHajjDMP <- subset(ElHajjDMP,ElHajjDMP$Gene%in%names(which(table(ElHajjDMP$Gene)>=2)))
CamaraschiDMP <- subset(CamaraschiDMP,CamaraschiDMP$Gene%in%names(which(table(CamaraschiDMP$Gene)>=2)))

ListDMP=list("YeungDMP"=YeungDMP,"BarberetDMP"=BarberetDMP,"DucreuxDMP"=DucreuxDMP,
             "TobiDMP"=TobiDMP,"NovakovicneoDMP"=NovakovicneoDMP,"MelamedDMP"=MelamedDMP,"KatariDMPcord"=KatariDMPcord,"KatariDMPplacenta"=KatariDMPplacenta,"ElHajjDMP"=ElHajjDMP,"CamaraschiDMP"=CamaraschiDMP)


for (k in 1:length(ListDMP)){
  cat(c(names(ListDMP[k])))
  cat(sep = '\n')
  cat(intersect(ListDMP[[k]]$Gene,IGlist$V1))
  cat(sep = '\n')
  cat("dim=")
  cat(length((intersect(ListDMP[[k]]$Gene,IGlist$V1))),sep = '\n')
  }

for (k in 1:length(ListDMP)){
  for (j in k:length(ListDMP)){
    if (names(ListDMP[k])!=names(ListDMP[j])){
      cat(c(names(ListDMP[k]),names(ListDMP[j])))
      cat(intersect(ListDMP[[k]]$CpG,ListDMP[[j]]$CpG),sep = '\n')
      cat(sep = '\n')
    }
  }
}  

for (k in 1:length(ListDMP)){
  for (j in k:length(ListDMP)){
    if (names(ListDMP[k])!=names(ListDMP[j])){
      cat(c(names(ListDMP[k]),names(ListDMP[j])))
      cat(sep = '\n')
      cat(intersect(ListDMP[[k]]$Gene,ListDMP[[j]]$Gene),sep = '\n')
      cat(sep = '\n')
    }
  }
}  

x=list("Yeung et al. (2021) neonates"=na.omit(unique(YeungDMP$Gene)),"Barberet et al. (2021)"=na.omit(unique(BarberetDMP$Gene)),
       "Tobi et al. (2021)"=na.omit(unique(TobiDMP$Gene)),"Novakovic et al. (2019) neonates"=na.omit(unique(NovakovicneoDMP$Gene)),
       "Melamed et al. (2015)"=na.omit(unique(MelamedDMP$Gene)),"Katari et al. (2009) cord blood"=na.omit(unique(KatariDMPcord$Gene)),"Katari et al. (2009) placenta"=na.omit(unique(KatariDMPplacenta$Gene)),
       "El Hajj et al. (2017)"=na.omit(unique(ElHajjDMP$Gene)),"Camaraschi et al. (2021)"=na.omit(unique(CamaraschiDMP$Gene)))

x=fromList(x)
names=c("Barberet et al. (2021)","Yeung et al. (2021) neonates","Novakovic et al. (2019) neonates","Camaraschi et al. (2021)","Tobi et al. (2021)","El Hajj et al. (2017)","Melamed et al. (2015)","Katari et al. (2009) cord blood","Katari et al. (2009) placenta")

x2=x[,names]
x_metadata = data.frame(
  set=c("Camaraschi et al. (2021)","Tobi et al. (2021)","Melamed et al. (2015)","Katari et al. (2009) cord blood","El Hajj et al. (2017)","Katari et al. (2009) placenta","Novakovic et al. (2019) neonates","Yeung et al. (2021) neonates","Barberet et al. (2021)"),
  celltype=c('Cord blood','Cord blood','Cord blood', 'Cord blood','Cord blood','Placenta','Bloodspots','Bloodspots','Buccal smears')
)

upset(x2,
      colnames(x2),name="",
      
      base_annotations=list('Intersection size'=intersection_size(bar_number_threshold = 1.2,text_colors =c(on_background = "black", on_bar = "black"),text = list(position = position_stack(vjust = 1)))+ ylim(c(0, 700))+ylab("Gene sets intersection size")),
      width_ratio=0.2,
      height_ratio=1,
      stripes=upset_stripes(geom=geom_segment(size=12),mapping=aes(color=celltype),colors=c('Buccal smears'='#4AACC6','Cord blood'='#92CDDC','Placenta'='#92CDDC','Bloodspots'='#92CDDC'),data=x_metadata),
      
      encode_sets=FALSE,
      matrix=(intersection_matrix(segment=geom_segment(linetype="solid",lineend="round",size=1.5,position = "identity"),outline_color=list(active=NA,inactive=NA))+ scale_color_manual(
        values=c('TRUE'='orange', 'FALSE'='white','Buccal smears'='#67d6c9','Cord blood'='#f1daa1','Placenta'='#f8c69d','Bloodspots'='#f0a592'),breaks=c('Placenta','Cord blood','Bloodspots','Buccal smears')
      )),
      themes=upset_default_themes(element_blank()),sort_sets='FALSE',
      set_sizes=(upset_set_size())+ geom_text(aes(label=..count..), hjust=1.4, stat='count')+ expand_limits(y=900)+ylab("Gene sets size (DMP)"),
      
      queries=list(
        upset_query(intersect='El Hajj et al. (2017)', color='azure3',fill='azure3'),
        upset_query(intersect='Katari et al. (2009) cord blood', color='azure3',fill='azure3'),
        upset_query(intersect='Katari et al. (2009) placenta', color='azure3',fill='azure3'),        
        upset_query(intersect="Barberet et al. (2021)", color='azure4',fill='azure4'),
        upset_query(intersect='Melamed et al. (2015)', color='azure3',fill='azure3'),
        
        upset_query(intersect=c("Barberet et al. (2021)","El Hajj et al. (2017)"), color='azure4',fill='azure4'),        
        upset_query(intersect=c('Katari et al. (2009) cord blood','Katari et al. (2009) placenta'), color='azure3',fill='azure3'),
        
        upset_query(intersect=c("Barberet et al. (2021)","Katari et al. (2009) cord blood"), color='azure4',fill='azure4'),
        upset_query(intersect='Novakovic et al. (2019) neonates', color='azure3',fill='azure3'),
        upset_query(intersect=c('Katari et al. (2009) placenta',"El Hajj et al. (2017)"), color='azure3',fill='azure3'),
        upset_query(intersect=c("Yeung et al. (2021) neonates","El Hajj et al. (2017)"), color='azure3',fill='azure3'),    
        
        upset_query(intersect=c('Yeung et al. (2021) neonates'), color='azure3',fill='azure3'),
        upset_query(intersect=c('Tobi et al. (2021)'), color='azure3',fill='azure3'),
        upset_query(intersect=c("Novakovic et al. (2019) neonates","El Hajj et al. (2017)"), color='azure3',fill='azure3'),  
        upset_query(intersect=c('Camaraschi et al. (2021)'), color='azure3',fill='azure3'),
        upset_query(intersect=c("Katari et al. (2009) cord blood","El Hajj et al. (2017)","Barberet et al. (2021)"), color='azure4',fill='azure4'),
        
        upset_query(intersect=c("Barberet et al. (2021)","Katari et al. (2009) cord blood","Melamed et al. (2015)"), color='azure4',fill='azure4'),
        upset_query(intersect=c("Barberet et al. (2021)","Katari et al. (2009) cord blood","Katari et al. (2009) placenta"), color='azure4',fill='azure4'),
        upset_query(intersect=c("Katari et al. (2009) cord blood","Katari et al. (2009) placenta","Barberet et al. (2021)","El Hajj et al. (2017)"), color='azure4',fill='azure4')
        
        ),
      
      intersections=list(c("Katari et al. (2009) placenta"),
                         c("Katari et al. (2009) cord blood"),c("Melamed et al. (2015)"), c("El Hajj et al. (2017)"),c("Tobi et al. (2021)"),c("Camaraschi et al. (2021)"),
                         c('Katari et al. (2009) cord blood','Katari et al. (2009) placenta'),c('Katari et al. (2009) placenta',"El Hajj et al. (2017)"),
                         c("Novakovic et al. (2019) neonates"),c("Yeung et al. (2021) neonates"),
                         c("Novakovic et al. (2019) neonates","El Hajj et al. (2017)"),c("Yeung et al. (2021) neonates","El Hajj et al. (2017)"),
                         c("Barberet et al. (2021)"),
                         c("Barberet et al. (2021)","El Hajj et al. (2017)"),
                         c("Barberet et al. (2021)","Katari et al. (2009) cord blood"),
                         c("Barberet et al. (2021)","El Hajj et al. (2017)","Katari et al. (2009) cord blood"),
                         c("Barberet et al. (2021)","Katari et al. (2009) cord blood","Melamed et al. (2015)"),
                         
                         c("Barberet et al. (2021)","Katari et al. (2009) cord blood","Katari et al. (2009) placenta"),
                         c("Katari et al. (2009) cord blood","Katari et al. (2009) placenta","Barberet et al. (2021)","El Hajj et al. (2017)")
                         )
      ,sort_intersections=FALSE)

#DMR analysis

ChenIVFDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Chen DMR IVF-ET hg38"))
ChenICSIDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Chen DMR ICSI-ET hg38"))
ChenIVFDMR = GRanges(seqnames=ChenIVFDMR$chr,ranges=IRanges(start=ChenIVFDMR$start,end=ChenIVFDMR$end),Deltabeta=ChenIVFDMR$Deltabeta)
ChenICSIDMR = GRanges(seqnames=ChenICSIDMR$chr,ranges=IRanges(start=ChenICSIDMR$start,end=ChenICSIDMR$end),Deltabeta=ChenICSIDMR$Deltabeta)
commonChenDMR=findOverlaps(ChenIVFDMR,ChenICSIDMR)

commonChenDMR=cbind(data.frame(ChenIVFDMR[commonChenDMR@from])[,c(1:3,6)],data.frame(ChenICSIDMR[commonChenDMR@to])[,6])
ChenDMR=subset(commonChenDMR,commonChenDMR[,4]==commonChenDMR[,5])
ChenDMR=ChenDMR[,-5]
ChenDMR = GRanges(seqnames=ChenDMR$seqnames,ranges=IRanges(start=ChenDMR$start,end=ChenDMR$end),Deltabeta=ChenDMR$Deltabeta)

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = rtracklayer::import.chain(path)

seqlevelsStyle(ChenDMR) = "UCSC"
ChenDMR = liftOver(ChenDMR, ch)
ChenDMR = unlist(ChenDMR)
genome(ChenDMR) = "hg19"
ChenDMR = new("gwaswloc", ChenDMR)

ChenDMR=data.frame(ChenDMR)
YeungDMRneo <- data.frame(read_excel("DMR annot.xlsx",sheet = "Yeung DMR IG neo"))
YeungDMRchild <- data.frame(read_excel("DMR annot.xlsx",sheet = "Yeung DMR IG child"))
BarberetDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Barberet DMR"))
DucreuxDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Ducreux DMR"))
NovakovicneoDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Novakovic DMR neonates"))
NovakovicaduDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Novakovic DMR adults"))
CastilloDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Castillo Fernandez DMR25%"))
EstillDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Estill DMR"))
CamprubiDMR <- data.frame(read_excel("DMR annot.xlsx",sheet = "Camprubi DMR"))

YeungDMRneo = GRanges(seqnames=YeungDMRneo$chr,ranges=IRanges(start=YeungDMRneo$start,end=YeungDMRneo$end))
YeungDMRchild = GRanges(seqnames=YeungDMRchild$chr,ranges=IRanges(start=YeungDMRchild$start,end=YeungDMRchild$end))
BarberetDMR = GRanges(seqnames=BarberetDMR$chr,ranges=IRanges(start=BarberetDMR$start,end=BarberetDMR$end))
DucreuxDMR = GRanges(seqnames=DucreuxDMR$chr,ranges=IRanges(start=DucreuxDMR$start,end=DucreuxDMR$end))
NovakovicneoDMR = GRanges(seqnames=NovakovicneoDMR$chr,ranges=IRanges(start=NovakovicneoDMR$start,end=NovakovicneoDMR$end))
NovakovicaduDMR = GRanges(seqnames=NovakovicaduDMR$chr,ranges=IRanges(start=NovakovicaduDMR$start,end=NovakovicaduDMR$end))
CastilloDMR = GRanges(seqnames=CastilloDMR$chr,ranges=IRanges(start=CastilloDMR$start,end=CastilloDMR$end))
EstillDMR = GRanges(seqnames=EstillDMR$chr,ranges=IRanges(start=EstillDMR$start,end=EstillDMR$end))
ChenDMR = GRanges(seqnames=ChenDMR$seqnames,ranges=IRanges(start=ChenDMR$start,end=ChenDMR$end),Deltabeta=ChenDMR$Deltabeta)
CamprubiDMR = GRanges(seqnames=CamprubiDMR$chr,ranges=IRanges(start=CamprubiDMR$start,end=CamprubiDMR$end),Deltabeta=CamprubiDMR$Deltabeta)


ListDMR=list("YeungDMRneo"=YeungDMRneo,"YeungDMRchild"=YeungDMRchild,"BarberetDMR"=BarberetDMR,"DucreuxDMR"=DucreuxDMR,
             "NovakovicneoDMR"=NovakovicneoDMR,"NovakovicaduDMR"=NovakovicaduDMR,"CastilloDMR"=CastilloDMR,
             "EstillDMR"=EstillDMR,"ChenDMR"=ChenDMR,"CamprubiDMR"=CamprubiDMR)

for (k in 1:length(ListDMR)){
  for (j in k:length(ListDMR)){
    if (names(ListDMR[k])!=names(ListDMR[j])){
      cat(c(names(ListDMR[k]),names(ListDMR[j])))
      cat(sep = '\n')
      a=findOverlaps(ListDMR[[k]], ListDMR[[j]],ignore.strand=TRUE)
      cat(sep = '\n')
      print(unique(ListDMR[[k]][a@from]))
      cat(sep = '\n')
      cat("dim=")
      cat(length(unique(findOverlaps(ListDMR[[k]], ListDMR[[j]],select="first",ignore.strand=TRUE))),sep = '\n')
    }
  }
}

IG_regions_start_end <- read_excel("IG regions start end.xlsx")
IG_regions_start_end=GRanges(seqnames=IG_regions_start_end$chr,ranges=IRanges(start=IG_regions_start_end$start,end=IG_regions_start_end$end),Region=IG_regions_start_end$Region)

for (k in 1:length(ListDMR)){
  cat(c(names(ListDMR[k])))
  cat(sep = '\n')
  a=findOverlaps(ListDMR[[k]],IG_regions_start_end,select="all",ignore.strand=TRUE)
  print(ListDMR[[k]][a@from])
  print(IG_regions_start_end[a@to]$Region)
  cat(sep = '\n')
} 

EstillDMR=matchGenes(EstillDMR,genes)
EstillDMR=data.frame(EstillDMR$name)
colnames(EstillDMR)="Gene"
EstillDMRGene <- data.frame(read_excel("DMR annot.xlsx",sheet = "Estill DMR IG genes"))
EstillDMR=data.frame(unique(c(EstillDMRGene$Gene,EstillDMR$Gene)))
colnames(EstillDMR)="Gene"

genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

YeungDMRneo$Gene=matchGenes(YeungDMRneo,genes)$name
YeungDMRchild$Gene=matchGenes(YeungDMRchild,genes)$name
BarberetDMR$Gene=matchGenes(BarberetDMR,genes)$name
DucreuxDMR$Gene=matchGenes(DucreuxDMR,genes)$name
NovakovicneoDMR$Gene=matchGenes(NovakovicneoDMR,genes)$name
NovakovicaduDMR$Gene=matchGenes(NovakovicaduDMR,genes)$name
CastilloDMR$Gene=matchGenes(CastilloDMR,genes)$name
CamprubiDMR$Gene=matchGenes(CamprubiDMR,genes)$name
ChenDMR$Gene=matchGenes(ChenDMR,genes)$name

ListDMRGene=list("YeungDMRneo"=YeungDMRneo,"YeungDMRchild"=YeungDMRchild,"BarberetDMR"=BarberetDMR,
                 "NovakovicneoDMR"=NovakovicneoDMR,"NovakovicaduDMR"=NovakovicaduDMR,"CastilloDMR"=CastilloDMR,
                 "EstillDMR"=EstillDMR,"ChenDMR"=ChenDMR,"Camprubi"=CamprubiDMR)

for (k in 1:length(ListDMRGene)){
  cat(c(names(ListDMRGene[k])))
  cat(sep = '\n')
  cat(intersect(ListDMRGene[[k]]$Gene,IGlist$V1))
  cat(sep = '\n')
} 

for (k in 1:length(ListDMRGene)){
  for (j in k:length(ListDMRGene)){
    if (names(ListDMRGene[k])!=names(ListDMRGene[j])){
      cat(c(names(ListDMRGene[k]),names(ListDMRGene[j])))
      cat(sep = '\n')
      cat(intersect(ListDMRGene[[k]]$Gene,ListDMRGene[[j]]$Gene),sep = '\n')
      cat(sep = '\n')
      cat("dim=")
      cat(length(intersect(ListDMRGene[[k]]$Gene,ListDMRGene[[j]]$Gene)),sep = '\n')
    }
  }
} 


x=list("Yeung et al. (2021) neonates"=na.omit(unique(YeungDMRneo$Gene)),"Yeung et al. (2021) children"=na.omit(unique(YeungDMRchild$Gene)),"Barberet et al. (2021)"=na.omit(unique(BarberetDMR$Gene)),
       "Novakovic et al. (2019) neonates"=na.omit(unique(NovakovicneoDMR$Gene)),"Novakovic et al. (2019) adults"=na.omit(unique(NovakovicaduDMR$Gene)),"Castillo-Fernandez et al. (2017)"=na.omit(unique(CastilloDMR$Gene)),
       "Estill et al. (2016)"=na.omit(unique(EstillDMR$Gene)),"Chen et al. (2020)"=na.omit(unique(ChenDMR$Gene)),"Camprubi et al. (2013)"=na.omit(unique(CamprubiDMR$Gene)))
x=fromList(x)
names=c("Novakovic et al. (2019) adults","Yeung et al. (2021) children","Barberet et al. (2021)","Yeung et al. (2021) neonates","Novakovic et al. (2019) neonates","Estill et al. (2016)","Chen et al. (2020)","Castillo-Fernandez et al. (2017)","Camprubi et al. (2013)")
x2=x[,names]
x_metadata = data.frame(
  set=c("Barberet et al. (2021)","Yeung et al. (2021) children","Yeung et al. (2021) neonates","Novakovic et al. (2019) adults","Novakovic et al. (2019) neonates","Castillo-Fernandez et al. (2017)","Estill et al. (2016)","Chen et al. (2020)","Camprubi et al. (2013)"),
  celltype=c('Buccal smears', 'Blood','Bloodspots','Blood','Bloodspots','Cord blood','Bloodspots','Cord blood','Placenta')
)


upset(x2,
      colnames(x2),name="",
      
      base_annotations=list('Intersection size'=intersection_size(bar_number_threshold = 0.67,text_colors =c(on_background = "black", on_bar = "black"),text = list(position = position_stack(vjust = 1.5)))+ ylim(c(0, 700))+ylab("Gene sets intersection size")),
      width_ratio=0.2,
      height_ratio=1,
      stripes=upset_stripes(geom=geom_segment(size=12),mapping=aes(color=celltype),colors=c('Buccal smears'='#4AACC6','Cord blood'='#92CDDC','Placenta'='#92CDDC','Bloodspots'='#92CDDC','Blood'='#30839A'),data=x_metadata),
      
      encode_sets=FALSE,
      matrix=(intersection_matrix(segment=geom_segment(linetype="solid",lineend="round",size=1.5,position = "identity"),outline_color=list(active=NA,inactive=NA))+ scale_color_manual(
        values=c('TRUE'='orange', 'FALSE'='white','Buccal smears'='#67d6c9','Cord blood'='#f1daa1','Same tissue intersection'='#f8c69d','Bloodspots'='#f0a592','Blood'="#81b2c5"),breaks=c('Same tissue intersection','Cord blood','Bloodspots','Buccal smears','Blood')
      )),
      themes=upset_default_themes(element_blank()),sort_sets='FALSE',
      set_sizes=(upset_set_size())+ geom_text(aes(label=..count..), hjust=1.4, stat='count')+ expand_limits(y=800)+ylab("Gene sets size (DMR)"),
      
      queries=list(
        upset_query(intersect='Chen et al. (2020)', color='azure3',fill='azure3'),
        upset_query(intersect='Estill et al. (2016)', color='azure3',fill='azure3'),
        upset_query(intersect='Novakovic et al. (2019) neonates', color='azure3',fill='azure3'),
        upset_query(intersect=c("Barberet et al. (2021)"), color='azure4',fill='azure4'),
        upset_query(intersect=c("Barberet et al. (2021)","Estill et al. (2016)"), color='azure4',fill='azure4'),
        upset_query(intersect=c("Estill et al. (2016)","Chen et al. (2020)"), color='azure3',fill='azure3'),
        upset_query(intersect=c('Yeung et al. (2021) neonates'), color='azure3',fill='azure3'),
        upset_query(intersect=c('Novakovic et al. (2019) adults',"Novakovic et al. (2019) neonates"), color='black',fill='black'),
        upset_query(intersect=c("Barberet et al. (2021)","Chen et al. (2020)"), color='azure4',fill='azure4'),
        upset_query(intersect=c('Estill et al. (2016)','Yeung et al. (2021) neonates'), color='azure3',fill='azure3'),
        upset_query(intersect='Novakovic et al. (2019) adults', color='black',fill='black'),
        upset_query(intersect=c('Chen et al. (2020)','Novakovic et al. (2019) neonates'), color='azure3',fill='azure3'),
        upset_query(intersect=c('Chen et al. (2020)','Estill et al. (2016)',"Barberet et al. (2021)","Camprubi et al. (2013)"), color='azure4',fill='azure4'),
        upset_query(intersect=c('Yeung et al. (2021) neonates','Chen et al. (2020)'), color='azure3',fill='azure3'),
        upset_query(intersect='Castillo-Fernandez et al. (2017)', color='azure3',fill='azure3'),
        upset_query(intersect=c('Camprubi et al. (2013)',"Estill et al. (2016)"), color='azure3',fill='azure3'),
        upset_query(intersect=c('Camprubi et al. (2013)'), color='azure3',fill='azure3'),
     upset_query(intersect=c('Yeung et al. (2021) children',"Yeung et al. (2021) neonates"), color='azure4',fill='azure4'),
        upset_query(intersect=c("Barberet et al. (2021)","Yeung et al. (2021) neonates"), color='azure4',fill='azure4')),
      
      intersections=list(c("Camprubi et al. (2013)"),
                         c("Castillo-Fernandez et al. (2017)"),c("Chen et al. (2020)"),
                         c("Estill et al. (2016)"),c("Novakovic et al. (2019) neonates"),c("Yeung et al. (2021) neonates"),
                         
                         c("Yeung et al. (2021) neonates","Estill et al. (2016)"),c("Estill et al. (2016)","Chen et al. (2020)"),
                         c("Estill et al. (2016)","Camprubi et al. (2013)"),
                         c("Novakovic et al. (2019) neonates","Chen et al. (2020)"),
                         c("Yeung et al. (2021) neonates","Chen et al. (2020)"),
                         
                         
                         
                         
                         
                         c("Barberet et al. (2021)"),c("Barberet et al. (2021)","Yeung et al. (2021) neonates"),c("Barberet et al. (2021)","Estill et al. (2016)"),c("Barberet et al. (2021)","Chen et al. (2020)"),c("Barberet et al. (2021)","Estill et al. (2016)","Chen et al. (2020)","Camprubi et al. (2013)"),
                         c("Yeung et al. (2021) children","Yeung et al. (2021) neonates"),
                         c("Novakovic et al. (2019) adults"),c("Novakovic et al. (2019) adults","Novakovic et al. (2019) neonates")),sort_intersections=FALSE)


library(clusterProfiler)

enrichGO(unique(c(CamprubiDMR$Gene)),keyType = "SYMBOL",'org.Hs.eg.db', ont="BP", qvalueCutoff=0.05)
enrichGO(unique(c(CastilloDMR$Gene,ChenDMR$Gene)),keyType = "SYMBOL",'org.Hs.eg.db', ont="BP", qvalueCutoff=0.05)
enrichGO(unique(c(YeungDMRneo$Gene,NovakovicneoDMR$Gene,EstillDMR$Gene)),keyType = "SYMBOL",'org.Hs.eg.db', ont="BP", qvalueCutoff=0.05)
enrichGO(unique(c(YeungDMRchild$Gene,NovakovicaduDMR$Gene)),keyType = "SYMBOL",'org.Hs.eg.db', ont="BP", qvalueCutoff=0.05)
enrichGO(unique(c(BarberetDMR$Gene,DucreuxDMR$Gene)),keyType = "SYMBOL",'org.Hs.eg.db', ont="BP", qvalueCutoff=0.05)


