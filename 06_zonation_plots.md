05\_zonation\_plots
================

``` r
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(seqinr)
library(tidyr)
```

Load saved go terms
-------------------

Established in [05\_annotation\_files](05_annotation_files.md)

``` r
load(file = "./raw_data/annotation.rds")
```

Import data
-----------

Generated in matlab\_03\_zonation\_reconstruction.m

``` r
library(readr)
zonation<-read_delim("./raw_data/zonation_table.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
```

    ## Parsed with column specification:
    ## cols(
    ##   `Gene name` = col_character(),
    ##   Crypt_mean = col_double(),
    ##   V1_mean = col_double(),
    ##   V2_mean = col_double(),
    ##   V3_mean = col_double(),
    ##   V4_mean = col_double(),
    ##   V5_mean = col_double(),
    ##   V6_mean = col_double(),
    ##   Crypt_SE = col_double(),
    ##   V1_SE = col_double(),
    ##   V2_SE = col_double(),
    ##   V3_SE = col_double(),
    ##   V4_SE = col_double(),
    ##   V5_SE = col_double(),
    ##   V6_SE = col_double(),
    ##   pval = col_double(),
    ##   qval = col_double()
    ## )

``` r
colnames(zonation) <- c("gene", 1,2,3,4,5,6,7,"se_crypt","se_1","se_2","se_3","se_4","se_5","se_6","pval","qval")
zon<-gather(dplyr::select(zonation,1:8),key=position,value=mean,-c(gene))
colnames(zonation) <- c("gene","mean_crypt","mean_1","mean_2","mean_3","mean_4","mean_5","mean_6",1,2,3,4,5,6,7,"pval","qval")
zon_se<-gather(dplyr::select(zonation,1,9:15),key=position,value=se,-c(gene))

zones<-inner_join(zon,zon_se, by = c("position","gene"))
zones$position <- as.numeric(as.character(zones$position))
```

Define GO term plot
-------------------

``` r
plot_go<-function(go,cutoff){

go_selection<-filter(gene_ids,name_1006==go)$external_gene_name

#filter go selection by expression value 
filtered_go_selection <- zones%>%
  group_by(gene)%>%
  summarise(max=max(mean))%>%
  filter(max>cutoff)%>%
  filter(gene %in% go_selection)

#scale expression matrix by max
rownames(zonation)<-zonation$gene
A1 = zonation%>%
  dplyr::select(2:8)
scaled<-A1/apply(A1,1,max)
scaled$gene<-rownames(scaled)
colnames(scaled) <- c(1,2,3,4,5,6,7,"gene")
scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%filtered_go_selection$gene)
#summarize by position:
plot_table<-scaled_zones%>%
  group_by(position)%>%
  summarise(value=mean(mean),sem=sd(mean)/sqrt(nrow(filter(scaled_zones,position==1))))

print(filtered_go_selection$gene)

#plot
ggplot(plot_table, aes(x=position, y=value))+
  geom_rect(xmin = -Inf, xmax = 1.5,   ymin = -Inf, ymax = Inf, fill =  brewer.pal(9,"Blues")[3] )+
  geom_ribbon(group = 1,aes(ymin = value-sem, ymax = value+sem), fill = brewer.pal(9,"Blues")[6]) +
  geom_line(group = 1,color=brewer.pal(9,"Blues")[9],size=1.5)+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),
     #              labels=c("Crypt","Villus Zone 1","Villus Zone 2", "Villus Zone 3","Villus Zone 4","Villus Zone 5", "Villus Zone 6"))+
     labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),axis.title.x=element_blank(),axis.title.y=element_blank())+
 # ylab("Expression relative to max")+
  ggtitle(paste(go,", ",gene_ids[gene_ids$name_1006==go,][1,6],sep=""))
}
```

Expression plots
----------------

``` r
#cluster 1
translation<-plot_go("translation",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##   [1] "Eef1a1"     "Eef1b2"     "Eef1d"      "Eef1g"      "Eef2"      
    ##   [6] "Eif1"       "Eif1a"      "Eif1ax"     "Eif2s2"     "Eif3a"     
    ##  [11] "Eif3b"      "Eif3c"      "Eif3e"      "Eif3f"      "Eif3g"     
    ##  [16] "Eif3h"      "Eif3i"      "Eif3k"      "Eif3l"      "Eif3m"     
    ##  [21] "Eif4a1"     "Eif4a2"     "Eif4b"      "Eif4e2"     "Eif4g1"    
    ##  [26] "Eif4h"      "Eif5"       "Eif5a"      "Eif5b"      "Eif6"      
    ##  [31] "Etf1"       "Fau"        "Gm10020"    "Gm10036"    "Gm10076"   
    ##  [36] "Gm11808"    "Gm2000"     "Gm6133"     "Gm9493"     "Gm9843"    
    ##  [41] "Mrpl12"     "Mrpl14"     "Mrpl15"     "Mrpl18"     "Mrpl20"    
    ##  [46] "Mrpl21"     "Mrpl23"     "Mrpl24"     "Mrpl30"     "Mrpl33"    
    ##  [51] "Mrpl34"     "Mrpl35"     "Mrpl36"     "Mrpl51"     "Mrpl52"    
    ##  [56] "Mrpl55"     "Mrpl57"     "Mrps12"     "Mrps14"     "Mrps15"    
    ##  [61] "Mrps16"     "Mrps18a"    "Mrps21"     "Mrps24"     "Mrps5"     
    ##  [66] "Mrps7"      "Nars"       "Nhp2"       "Rbm3"       "Rpl10"     
    ##  [71] "Rpl10-ps3"  "Rpl10a"     "Rpl11"      "Rpl12"      "Rpl13"     
    ##  [76] "Rpl13a"     "Rpl13a-ps1" "Rpl14"      "Rpl15"      "Rpl17"     
    ##  [81] "Rpl18"      "Rpl18a"     "Rpl19"      "Rpl21"      "Rpl22"     
    ##  [86] "Rpl22l1"    "Rpl23"      "Rpl23a"     "Rpl23a-ps3" "Rpl24"     
    ##  [91] "Rpl26"      "Rpl27"      "Rpl27-ps3"  "Rpl27a"     "Rpl28"     
    ##  [96] "Rpl29"      "Rpl3"       "Rpl30"      "Rpl31"      "Rpl32"     
    ## [101] "Rpl34"      "Rpl35"      "Rpl35a"     "Rpl36"      "Rpl36a"    
    ## [106] "Rpl36al"    "Rpl37"      "Rpl37a"     "Rpl38"      "Rpl39"     
    ## [111] "Rpl4"       "Rpl41"      "Rpl5"       "Rpl6"       "Rpl7"      
    ## [116] "Rpl7a"      "Rpl8"       "Rpl9"       "Rpl9-ps6"   "Rps11"     
    ## [121] "Rps12"      "Rps12-ps3"  "Rps13"      "Rps14"      "Rps15"     
    ## [126] "Rps15a"     "Rps16"      "Rps17"      "Rps18"      "Rps19"     
    ## [131] "Rps2"       "Rps20"      "Rps21"      "Rps23"      "Rps24"     
    ## [136] "Rps26"      "Rps27"      "Rps27a"     "Rps27l"     "Rps27rt"   
    ## [141] "Rps28"      "Rps29"      "Rps3"       "Rps3a1"     "Rps4x"     
    ## [146] "Rps5"       "Rps6"       "Rps7"       "Rps8"       "Rps9"      
    ## [151] "Rpsa"       "Rrbp1"      "Uba52"

``` r
splicing<-plot_go("RNA splicing",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Alyref"    "Celf1"     "Ddx39b"    "Dhx15"     "Gm8186"   
    ##  [6] "Hnrnpa1"   "Hnrnpa2b1" "Hnrnpa3"   "Hnrnpc"    "Hnrnpf"   
    ## [11] "Hnrnph1"   "Hnrnpk"    "Hnrnpm"    "Hnrnpu"    "Hspa8"    
    ## [16] "Lsm2"      "Lsm3"      "Lsm5"      "Lsm6"      "Luc7l3"   
    ## [21] "Mbnl2"     "Nono"      "Pabpc1"    "Phf5a"     "Ptbp3"    
    ## [26] "Raly"      "Rbm25"     "Rbm39"     "Rbm8a"     "Rp9"      
    ## [31] "Sf3b5"     "Sf3b6"     "Sfpq"      "Snrpb"     "Snrpd1"   
    ## [36] "Snrpd2"    "Snrpd3"    "Snrpe"     "Snrpf"     "Snrpg"    
    ## [41] "Son"       "Srrm2"     "Srsf1"     "Srsf11"    "Srsf2"    
    ## [46] "Srsf3"     "Srsf5"     "Srsf6"     "Srsf7"     "Tardbp"   
    ## [51] "Thoc7"     "Tra2b"     "U2af1"     "Ybx1"

``` r
trasncription<-plot_go("transcription, DNA-templated",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##   [1] "2700060E02Rik" "Aes"           "Aff4"          "Ago2"         
    ##   [5] "Anp32a"        "Arid1a"        "Atf3"          "Atf4"         
    ##   [9] "Atxn1"         "Bach1"         "Bhlhe40"       "Birc5"        
    ##  [13] "Btf3"          "Btg2"          "Bzw1"          "C1d"          
    ##  [17] "Casz1"         "Ccnd1"         "Ccnl1"         "Cdc73"        
    ##  [21] "Cdx1"          "Cdx2"          "Cebpg"         "Cggbp1"       
    ##  [25] "Chchd2"        "Chchd3"        "Chmp1a"        "Churc1"       
    ##  [29] "Cidec"         "Cited2"        "Cnbp"          "Creb3l3"      
    ##  [33] "Csnk2a1"       "Ctbp1"         "Ctnnb1"        "Ctnnd1"       
    ##  [37] "Ddi2"          "Ddit3"         "Ddx3x"         "Dedd2"        
    ##  [41] "Dpy30"         "Drap1"         "Dynll1"        "Edf1"         
    ##  [45] "Elf3"          "Ell2"          "Elof1"         "Eny2"         
    ##  [49] "Esrra"         "Ets2"          "Ewsr1"         "Fosl2"        
    ##  [53] "Foxo3"         "Fubp1"         "Gata6"         "Gatad2a"      
    ##  [57] "Gtf2a2"        "Gtf2h5"        "Hdgf"          "Hif1a"        
    ##  [61] "Hint1"         "Hmg20b"        "Hmgb2"         "Hnf4a"        
    ##  [65] "Hnf4g"         "Hnrnpab"       "Hnrnpd"        "Hnrnpk"       
    ##  [69] "Hspa8"         "Id1"           "Id3"           "Ier5"         
    ##  [73] "Irf1"          "Irf2bp2"       "Irf7"          "Jmjd1c"       
    ##  [77] "Jun"           "Junb"          "Jund"          "Klf3"         
    ##  [81] "Klf4"          "Klf5"          "Litaf"         "Lpin2"        
    ##  [85] "Lrrfip1"       "Maf"           "Mafb"          "Mapk1"        
    ##  [89] "Mapk13"        "Mapk3"         "Med13"         "Med28"        
    ##  [93] "Mlxipl"        "Morf4l1"       "Morf4l2"       "Mxi1"         
    ##  [97] "Naca"          "Nfe2l2"        "Nfib"          "Nono"         
    ## [101] "Nr1i2"         "Nr3c1"         "Nrip1"         "Onecut2"      
    ## [105] "Pa2g4"         "Parp14"        "Pbrm1"         "Pcbd1"        
    ## [109] "Phb2"          "Phf5a"         "Pkn2"          "Polr1d"       
    ## [113] "Polr2e"        "Polr2f"        "Polr2j"        "Polr2k"       
    ## [117] "Polr2l"        "Ppp1r1b"       "Prr13"         "Pura"         
    ## [121] "Purb"          "Raly"          "Rbbp7"         "Rbm39"        
    ## [125] "Rnf10"         "Rps3"          "Sap30l"        "Sec14l2"      
    ## [129] "Sertad1"       "Sfpq"          "Slirp"         "Smad4"        
    ## [133] "Smad7"         "Snf8"          "Son"           "Sra1"         
    ## [137] "Srsf5"         "Ssrp1"         "Stat3"         "Sub1"         
    ## [141] "Taf10"         "Tardbp"        "Tbx3"          "Tcf25"        
    ## [145] "Tcf7l2"        "Tle3"          "Tsc22d1"       "Tsc22d4"      
    ## [149] "Txn1"          "Ube2l3"        "Vdr"           "Wasl"         
    ## [153] "Xbp1"          "Ybx1"          "Ybx3"          "Zbtb7a"       
    ## [157] "Zbtb7b"        "Zfp703"        "Zfpm1"         "Zmiz1"

``` r
#cluster 2
mito<-plot_go("mitochondrion",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##   [1] "1110008F13Rik" "2010107E04Rik" "2410015M20Rik" "AA467197"     
    ##   [5] "Abcd1"         "Abcd3"         "Acaa1a"        "Acaa1b"       
    ##   [9] "Acaa2"         "Acadl"         "Acadm"         "Acadvl"       
    ##  [13] "Acat1"         "Aco1"          "Aco2"          "Acot13"       
    ##  [17] "Acox1"         "Acsl3"         "Acsl5"         "Adh1"         
    ##  [21] "Ak2"           "Ak3"           "Akr1b7"        "Akr7a5"       
    ##  [25] "Aldh18a1"      "Aldh1b1"       "Aldh2"         "Aldh4a1"      
    ##  [29] "Ap2m1"         "Arg2"          "Arglu1"        "Asah2"        
    ##  [33] "Atp5a1"        "Atp5b"         "Atp5c1"        "Atp5d"        
    ##  [37] "Atp5e"         "Atp5f1"        "Atp5g1"        "Atp5g2"       
    ##  [41] "Atp5g3"        "Atp5h"         "Atp5j"         "Atp5j2"       
    ##  [45] "Atp5k"         "Atp5l"         "Atp5o"         "Atpif1"       
    ##  [49] "Aurkaip1"      "Bad"           "Bak1"          "Bcap31"       
    ##  [53] "Bnip3l"        "Bola3"         "Bri3bp"        "Bsg"          
    ##  [57] "C1qbp"         "Casp1"         "Casp8"         "Cat"          
    ##  [61] "Cct7"          "Chchd1"        "Chchd10"       "Chchd2"       
    ##  [65] "Chchd3"        "Chchd4"        "Chchd7"        "Cisd1"        
    ##  [69] "Cisd2"         "Cisd3"         "Ckb"           "Ckmt1"        
    ##  [73] "Clic1"         "Clic4"         "Cltc"          "Cmc2"         
    ##  [77] "Coa3"          "Coa5"          "Comt"          "Cox14"        
    ##  [81] "Cox17"         "Cox20"         "Cox4i1"        "Cox5a"        
    ##  [85] "Cox6a1"        "Cox6b1"        "Cox6b2"        "Cox6c"        
    ##  [89] "Cox7a1"        "Cox7a2"        "Cox7a2l"       "Cox7b"        
    ##  [93] "Cox7c"         "Cox8a"         "Cps1"          "Cpt2"         
    ##  [97] "Crat"          "Crot"          "Cs"            "Ctsb"         
    ## [101] "Cyb5a"         "Cyb5b"         "Cyb5r3"        "Cyba"         
    ## [105] "Cyc1"          "Cycs"          "Dbi"           "Ddx3x"        
    ## [109] "Dgat2"         "Dhrs1"         "Dhrs4"         "Dlst"         
    ## [113] "Dnaja1"        "Dnajc15"       "Dnajc19"       "Dnajc5"       
    ## [117] "Dtymk"         "Dut"           "Dynll1"        "Ech1"         
    ## [121] "Eny2"          "Etfa"          "Etfb"          "Fahd1"        
    ## [125] "Fam136a"       "Fam162a"       "Fam213a"       "Fis1"         
    ## [129] "Fkbp4"         "Fkbp8"         "Foxo3"         "Fth1"         
    ## [133] "Gadd45gip1"    "Gapdh"         "Ghitm"         "Gk"           
    ## [137] "Glrx"          "Glrx5"         "Gls"           "Glud1"        
    ## [141] "Gm10250"       "Gng5"          "Golph3"        "Gpd1"         
    ## [145] "Gpd2"          "Gpx1"          "Gpx4"          "Grpel1"       
    ## [149] "Gsk3b"         "Gsr"           "Gstk1"         "Gstp1"        
    ## [153] "Hadh"          "Hadha"         "Hadhb"         "Higd1a"       
    ## [157] "Higd2a"        "Hint2"         "Hjurp"         "Hkdc1"        
    ## [161] "Hsd17b10"      "Hsd17b4"       "Hsp90ab1"      "Hspa5"        
    ## [165] "Hspa9"         "Hspd1"         "Hspe1"         "Idh1"         
    ## [169] "Idh3a"         "Idh3b"         "Idh3g"         "Immt"         
    ## [173] "Jtb"           "Kras"          "Lap3"          "Ldha"         
    ## [177] "Letm1"         "Lipe"          "Lypla1"        "Maoa"         
    ## [181] "Maob"          "Map1lc3b"      "Map2k1"        "Map2k2"       
    ## [185] "Mapk1"         "Mapk3"         "Mapk8"         "Marc2"        
    ## [189] "March5"        "Mcl1"          "Mcu"           "Mdh1"         
    ## [193] "Mdh2"          "Me2"           "Mfn2"          "Mfsd7b"       
    ## [197] "Mgst1"         "Micu1"         "Minos1"        "Mpc1"         
    ## [201] "Mpc2"          "Mrpl12"        "Mrpl14"        "Mrpl15"       
    ## [205] "Mrpl18"        "Mrpl20"        "Mrpl21"        "Mrpl23"       
    ## [209] "Mrpl24"        "Mrpl28"        "Mrpl30"        "Mrpl33"       
    ## [213] "Mrpl34"        "Mrpl35"        "Mrpl36"        "Mrpl42"       
    ## [217] "Mrpl51"        "Mrpl52"        "Mrpl54"        "Mrpl55"       
    ## [221] "Mrpl57"        "Mrps12"        "Mrps14"        "Mrps15"       
    ## [225] "Mrps16"        "Mrps21"        "Mrps24"        "Mrps26"       
    ## [229] "Mrps33"        "Mrps36"        "Mrps5"         "Mrps7"        
    ## [233] "Msra"          "mt-Atp6"       "mt-Atp8"       "mt-Co1"       
    ## [237] "mt-Co2"        "mt-Co3"        "mt-Cytb"       "mt-Nd1"       
    ## [241] "mt-Nd2"        "mt-Nd3"        "mt-Nd4"        "mt-Nd4l"      
    ## [245] "mt-Nd5"        "Mtch1"         "Mtch2"         "Mxd1"         
    ## [249] "Nars"          "Ndfip2"        "Ndufa1"        "Ndufa10"      
    ## [253] "Ndufa12"       "Ndufa13"       "Ndufa2"        "Ndufa3"       
    ## [257] "Ndufa4"        "Ndufa5"        "Ndufa6"        "Ndufa7"       
    ## [261] "Ndufa8"        "Ndufa9"        "Ndufab1"       "Ndufaf4"      
    ## [265] "Ndufb10"       "Ndufb11"       "Ndufb2"        "Ndufb3"       
    ## [269] "Ndufb4"        "Ndufb5"        "Ndufb6"        "Ndufb7"       
    ## [273] "Ndufb8"        "Ndufb9"        "Ndufc1"        "Ndufc2"       
    ## [277] "Ndufs2"        "Ndufs3"        "Ndufs4"        "Ndufs5"       
    ## [281] "Ndufs6"        "Ndufs7"        "Ndufs8"        "Ndufv1"       
    ## [285] "Ndufv2"        "Ndufv3"        "Nipsnap3b"     "Nme1"         
    ## [289] "Nol7"          "Nt5c"          "Nudt19"        "Nxf1"         
    ## [293] "Oat"           "Ociad2"        "Ogdh"          "Otc"          
    ## [297] "Park7"         "Pccb"          "Pdha1"         "Pebp1"        
    ## [301] "Perp"          "Pet100"        "Phb"           "Phb2"         
    ## [305] "Phyh"          "Pisd"          "Pkm"           "Pmaip1"       
    ## [309] "Pnkd"          "Pnpla7"        "Pold4"         "Pon2"         
    ## [313] "Por"           "Ppp1cc"        "Prdx1"         "Prdx2"        
    ## [317] "Prdx3"         "Prdx5"         "Prelid1"       "Prkca"        
    ## [321] "Psap"          "Ptrh1"         "Pycard"        "Qdpr"         
    ## [325] "Rab11a"        "Rab24"         "Rfk"           "Rhoa"         
    ## [329] "Rmdn3"         "Rnasel"        "Rnf5"          "Romo1"        
    ## [333] "Rpl34"         "Rpl35a"        "Rps14"         "Rps15a"       
    ## [337] "Rps3"          "Sco2"          "Scp2"          "Sdha"         
    ## [341] "Sdhb"          "Sdhc"          "Sdhd"          "Sfxn1"        
    ## [345] "Sgk1"          "Sh3glb1"       "Slc25a10"      "Slc25a11"     
    ## [349] "Slc25a15"      "Slc25a22"      "Slc25a24"      "Slc25a3"      
    ## [353] "Slc25a36"      "Slc25a39"      "Slc25a45"      "Slc25a5"      
    ## [357] "Slc25a51"      "Slirp"         "Smdt1"         "Sod1"         
    ## [361] "Sod2"          "Sord"          "Src"           "Sri"          
    ## [365] "Stard7"        "Stat3"         "Sucla2"        "Suclg1"       
    ## [369] "Suclg2"        "Synj2bp"       "Timm10b"       "Timm13"       
    ## [373] "Timm8b"        "Tmem135"       "Tmem14c"       "Tmem160"      
    ## [377] "Tnfrsf1a"      "Tomm20"        "Tomm22"        "Tomm40"       
    ## [381] "Tomm5"         "Tomm7"         "Tomm70a"       "Trim31"       
    ## [385] "Tspo"          "Tstd1"         "Txn1"          "Txn2"         
    ## [389] "Txnrd1"        "Ubb"           "Ucp2"          "Uqcc2"        
    ## [393] "Uqcr10"        "Uqcr11"        "Uqcrb"         "Uqcrc1"       
    ## [397] "Uqcrc2"        "Uqcrfs1"       "Uqcrh"         "Uqcrq"        
    ## [401] "Usmg5"         "Vat1"          "Vdac1"         "Vdac2"        
    ## [405] "Vdac3"         "Ykt6"          "Ywhae"         "Ywhaz"

``` r
glutathione<-plot_go("glutathione transferase activity",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Eef1g" "Gsta1" "Gstk1" "Gstm1" "Gstm3" "Gsto1" "Gstp1" "Gstt1"
    ##  [9] "Ltc4s" "Mgst1" "Mgst2" "Mgst3"

``` r
acute<-plot_go("acute-phase response",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "Reg3a" "Reg3b" "Reg3g" "Stat3"

``` r
#cluster 3
absorp<-plot_go("intestinal absorption",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "Aco1"   "Cd36"   "Ddi2"   "F11r"   "Fabp1"  "Mogat2" "Vdr"

``` r
ion<-plot_go("ion transport",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Ano6"     "Atox1"    "Atp1a1"   "Atp1b1"   "Atp2a2"   "Atp2b1"  
    ##  [7] "Atp5a1"   "Atp5b"    "Atp5c1"   "Atp5d"    "Atp5e"    "Atp5f1"  
    ## [13] "Atp5g1"   "Atp5g2"   "Atp5g3"   "Atp5h"    "Atp5j"    "Atp5j2"  
    ## [19] "Atp5k"    "Atp5l"    "Atp5o"    "Atp6v0a2" "Atp6v0b"  "Atp6v0e" 
    ## [25] "Atp6v1f"  "Atp6v1g1" "Clcn3"    "Cldn15"   "Clic1"    "Clic4"   
    ## [31] "Clic5"    "Gm10250"  "Gpr89"    "Heph"     "Kcne3"    "Kcnk5"   
    ## [37] "Lasp1"    "Letm1"    "Lrrc8b"   "Mcu"      "Micu1"    "mt-Atp6" 
    ## [43] "mt-Atp8"  "Nipa2"    "Nipal1"   "Sfxn1"    "Slc12a2"  "Slc17a4" 
    ## [49] "Slc26a2"  "Slc31a1"  "Slc34a2"  "Slc39a4"  "Slc39a5"  "Slc39a9" 
    ## [55] "Slc4a4"   "Slc4a7"   "Slc5a12"  "Slc6a8"   "Slc9a2"   "Smdt1"   
    ## [61] "Tcn2"     "Tmco1"    "Tmem37"   "Tomm40"   "Tspo"     "Ttyh2"   
    ## [67] "Vdac1"    "Vdac2"    "Vdac3"

``` r
carboxy<-plot_go("carboxypeptidase activity",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "Ace"      "Ace2"     "Cndp2"    "Ctsz"     "Naaladl1"

``` r
#cluster 4
brush<-plot_go("brush border",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Actn4"   "Actr3"   "Anks4b"  "Aqp1"    "Capza2"  "Capzb"   "Cdhr2"  
    ##  [8] "Coro2a"  "Cubn"    "Dab1"    "Ddi2"    "Diaph1"  "Enpep"   "Eps8"   
    ## [15] "Espn"    "Ezr"     "Lima1"   "Mme"     "Myh14"   "Myl12b"  "Myl6"   
    ## [22] "Myo15b"  "Myo18a"  "Myo1a"   "Myo1d"   "Myo1e"   "Myo7b"   "Plec"   
    ## [29] "Pls1"    "Sis"     "Slc15a1" "Slc2a2"  "Slc34a2" "Slc9a2"  "Soat2"  
    ## [36] "Treh"    "Ush1c"   "Vil1"

``` r
plasma<-plot_go("plasma membrane",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##   [1] "1810011O10Rik" "1810055G02Rik" "2610528J11Rik" "AA986860"     
    ##   [5] "Abcc2"         "Abhd2"         "Abi1"          "Abr"          
    ##   [9] "Ace"           "Ace2"          "Acox1"         "Actb"         
    ##  [13] "Actg1"         "Ada"           "Adap1"         "Adh1"         
    ##  [17] "Adipor1"       "Adipor2"       "Adss"          "Agpat3"       
    ##  [21] "Ahcyl1"        "Ahnak"         "Alpi"          "Amn"          
    ##  [25] "Ankrd11"       "Anks4b"        "Ano6"          "Anpep"        
    ##  [29] "Anxa2"         "Anxa4"         "Ap2a2"         "Ap2m1"        
    ##  [33] "Ap2s1"         "Aplp1"         "App"           "Aqp1"         
    ##  [37] "Aqp7"          "Arf1"          "Arf4"          "Arf5"         
    ##  [41] "Arf6"          "Arhgef5"       "Arl4a"         "Arpc2"        
    ##  [45] "Asah2"         "Atf4"          "Atp10b"        "Atp1a1"       
    ##  [49] "Atp1b1"        "Atp2b1"        "Atp5a1"        "Atp5b"        
    ##  [53] "Atp5o"         "Atp6v0a2"      "Atp6v1g1"      "Atp8b1"       
    ##  [57] "AU040320"      "B2m"           "B4galnt1"      "Baiap2l1"     
    ##  [61] "Baiap2l2"      "Bcar1"         "Bsg"           "C1qbp"        
    ##  [65] "C2cd2l"        "Camk2d"        "Camk2n1"       "Car4"         
    ##  [69] "Casp8"         "Cat"           "Cct3"          "Cd164"        
    ##  [73] "Cd24a"         "Cd2ap"         "Cd36"          "Cd47"         
    ##  [77] "Cd55"          "Cd63"          "Cd74"          "Cd81"         
    ##  [81] "Cdc42"         "Cdc42bpb"      "Cdc42ep5"      "Cdc42se2"     
    ##  [85] "Cdh1"          "Cdh17"         "Cdhr2"         "Cdhr5"        
    ##  [89] "Cdipt"         "Ceacam1"       "Ceacam20"      "Cfl1"         
    ##  [93] "Cgn"           "Chmp2b"        "Chmp3"         "Chp1"         
    ##  [97] "Cib1"          "Clcn3"         "Cldn15"        "Cldn23"       
    ## [101] "Cldn3"         "Cldn4"         "Cldn7"         "Clec2d"       
    ## [105] "Clec2h"        "Clic1"         "Clic4"         "Cltb"         
    ## [109] "Cobl"          "Comt"          "Copb1"         "Coro1b"       
    ## [113] "Coro1c"        "Crb3"          "Csde1"         "Csnk2b"       
    ## [117] "Ctnna1"        "Ctnnb1"        "Ctnnd1"        "Cttn"         
    ## [121] "Cubn"          "Cxadr"         "Cyba"          "Dag1"         
    ## [125] "Ddi2"          "Diaph1"        "Dnajc5"        "Dnm2"         
    ## [129] "Dpep1"         "Dpp4"          "Dram2"         "Dsc2"         
    ## [133] "Dsg2"          "Dsp"           "Eef1a1"        "Eef2"         
    ## [137] "Efna1"         "Efr3a"         "Ehd1"          "Eif5"         
    ## [141] "Emp1"          "Eno1"          "Enpep"         "Entpd8"       
    ## [145] "Epb41l3"       "Epb41l4b"      "Epcam"         "Epn1"         
    ## [149] "Eps8"          "Eps8l2"        "Errfi1"        "Esyt2"        
    ## [153] "Etnk1"         "Ets2"          "Ewsr1"         "Ezr"          
    ## [157] "F11r"          "Fam120a"       "Fat1"          "Fbp2"         
    ## [161] "Frk"           "G3bp1"         "Gabarap"       "Gapdh"        
    ## [165] "Gde1"          "Gdpd2"         "Ggt1"          "Glo1"         
    ## [169] "Gna11"         "Gna13"         "Gnai2"         "Gnai3"        
    ## [173] "Gnas"          "Gnb1"          "Gng12"         "Gng5"         
    ## [177] "Golph3"        "Gpi1"          "Gprc5a"        "Gsdmd"        
    ## [181] "Gsk3b"         "Gstp1"         "Gucy2c"        "H2-Aa"        
    ## [185] "H2-Ab1"        "H2-D1"         "H2-Q1"         "H2-Q2"        
    ## [189] "H2-T23"        "Hbegf"         "Heph"          "Hint1"        
    ## [193] "Hmgb1"         "Hras"          "Hsp90aa1"      "Hsp90ab1"     
    ## [197] "Hsp90b1"       "Hspa5"         "Hspa8"         "Hspd1"        
    ## [201] "Ifitm2"        "Ifitm3"        "Ifngr1"        "Ifngr2"       
    ## [205] "Igsf9"         "Il17rc"        "Iqgap1"        "Itga6"        
    ## [209] "Itm2b"         "Jak2"          "Jup"           "Kcne3"        
    ## [213] "Kras"          "Krt19"         "Lamp1"         "Lamp2"        
    ## [217] "Lamtor1"       "Lipe"          "Litaf"         "Lrrc8b"       
    ## [221] "Lrrfip1"       "Lsr"           "Ly6g6c"        "Lypd8"        
    ## [225] "Mal2"          "Mall"          "Map2k1"        "Mapk1"        
    ## [229] "Mapk3"         "Mark2"         "Mdm2"          "Mep1b"        
    ## [233] "Metap2"        "Mfsd7b"        "Mgst2"         "Mien1"        
    ## [237] "Misp"          "Mme"           "Morf4l2"       "Mpp5"         
    ## [241] "Mrpl42"        "Muc13"         "Naaladl1"      "Nckap1"       
    ## [245] "Net1"          "Nfe2l2"        "Nfkbia"        "Nipa2"        
    ## [249] "Nlrp6"         "Npc1"          "Nt5e"          "Ocln"         
    ## [253] "Ogt"           "P4hb"          "Park7"         "Pcdh1"        
    ## [257] "Pcyt1a"        "Pdcd10"        "Perp"          "Phb"          
    ## [261] "Pi4k2b"        "Pigr"          "Pkp3"          "Plek2"        
    ## [265] "Plekha1"       "Plp2"          "Plpp1"         "Plxnb2"       
    ## [269] "Pmp22"         "Pon2"          "Ppm1a"         "Ppp1ca"       
    ## [273] "Ppp2ca"        "Prkar2a"       "Prkca"         "Prkcd"        
    ## [277] "Prkcz"         "Prss32"        "Psenen"        "Ptger4"       
    ## [281] "Ptp4a1"        "Ptp4a2"        "Ptprh"         "Ptprr"        
    ## [285] "Pttg1ip"       "Rab10"         "Rab11a"        "Rab18"        
    ## [289] "Rab25"         "Rab3gap2"      "Rab5c"         "Rab8a"        
    ## [293] "Rac1"          "Rala"          "Ralgps2"       "Rap1a"        
    ## [297] "Rap2c"         "Raph1"         "Rer1"          "Rhoa"         
    ## [301] "Rhob"          "Rhoc"          "Rhod"          "Rhof"         
    ## [305] "Rhou"          "Rnpep"         "Rps3"          "Rpsa"         
    ## [309] "S100a6"        "Samhd1"        "Scamp5"        "Sdcbp"        
    ## [313] "Sdhb"          "Sema4g"        "Sema6a"        "Serbp1"       
    ## [317] "Serinc3"       "Sgk1"          "Sh3d19"        "Slc12a2"      
    ## [321] "Slc16a10"      "Slc16a3"       "Slc17a4"       "Slc26a2"      
    ## [325] "Slc26a3"       "Slc27a4"       "Slc2a2"        "Slc31a1"      
    ## [329] "Slc34a2"       "Slc35g1"       "Slc39a4"       "Slc39a5"      
    ## [333] "Slc3a1"        "Slc3a2"        "Slc43a2"       "Slc46a1"      
    ## [337] "Slc4a4"        "Slc4a7"        "Slc51a"        "Slc51b"       
    ## [341] "Slc52a2"       "Slc5a12"       "Slc6a19"       "Slc6a20a"     
    ## [345] "Slc6a6"        "Slc7a7"        "Slc7a8"        "Slc7a9"       
    ## [349] "Slc9a2"        "Slc9a3r1"      "Smad7"         "Smagp"        
    ## [353] "Smap1"         "Smpdl3b"       "Snf8"          "Sod1"         
    ## [357] "Spint1"        "Spint2"        "Sppl2a"        "Sptbn1"       
    ## [361] "Sra1"          "Src"           "Sri"           "Ssfa2"        
    ## [365] "St14"          "Stap2"         "Stat3"         "Stk17b"       
    ## [369] "Stx3"          "Stx7"          "Suclg1"        "Suclg2"       
    ## [373] "Sumo1"         "Susd2"         "Svip"          "Synpo"        
    ## [377] "Tax1bp3"       "Tfrc"          "Tgoln1"        "Tm4sf20"      
    ## [381] "Tmed10"        "Tmem184a"      "Tmem30b"       "Tmem59"       
    ## [385] "Tnfrsf1a"      "Tonsl"         "Tradd"         "Treh"         
    ## [389] "Tspan13"       "Tspan15"       "Ttyh2"         "Ube2c"        
    ## [393] "Ube2d3"        "Ubl3"          "Ubqln1"        "Vamp8"        
    ## [397] "Vapa"          "Vasp"          "Vdac1"         "Vil1"         
    ## [401] "Vmp1"          "Vnn1"          "Vps37b"        "Xpnpep2"      
    ## [405] "Ywhae"         "Ywhah"

``` r
proj<-plot_go("cell projection organization",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "Actr3"      "Asap1"      "Capzb"      "Cfl1"       "D1Ertd622e"
    ## [6] "Ehd1"       "Ift20"      "Myo1a"      "Rab8a"

``` r
#cluster 5
adhesion<-plot_go("cell adhesion",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Ada"     "App"     "Atp1b1"  "Bcar1"   "Cd164"   "Cd24a"   "Cd36"   
    ##  [8] "Cd47"    "Cd63"    "Cdh1"    "Cdh17"   "Cdhr2"   "Cdhr5"   "Cib1"   
    ## [15] "Ctnna1"  "Ctnnb1"  "Ctnnd1"  "Cxadr"   "Dpp4"    "Dsc2"    "Dsg2"   
    ## [22] "F11r"    "Fat1"    "Itga6"   "Jup"     "Lama3"   "Lamb3"   "Lpp"    
    ## [29] "Pcdh1"   "Perp"    "Pkn2"    "Pkp3"    "Ptprf"   "Rac1"    "Rhoa"   
    ## [36] "Rhob"    "Specc1l" "Src"     "Tgfbi"   "Vmp1"

``` r
junction<-plot_go("cell-cell junction",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ##  [1] "Actn4"    "Actr3"    "Ahnak"    "App"      "Aqp7"     "Cd2ap"   
    ##  [7] "Cdc42"    "Cdc42bpb" "Ceacam1"  "Cfl1"     "Clic4"    "Ctnna1"  
    ## [13] "Ctnnb1"   "Ctnnd1"   "Cxadr"    "Dag1"     "Dsp"      "Epb41l3" 
    ## [19] "F11r"     "Fat1"     "Iqgap1"   "Jup"      "Krt8"     "Map2k2"  
    ## [25] "Myo1e"    "Ocln"     "Pcdh1"    "Prkcd"    "Prkcz"    "Slc2a2"  
    ## [31] "Slc5a1"   "Stx3"     "Twf1"     "Wdr1"

``` r
lipo<-plot_go("lipoprotein biosynthetic process",1E-4)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "Apoa1"   "Apob"    "Apobec1"

``` r
clusters_plot<-plot_grid(translation, splicing, trasncription,mito, glutathione, acute,absorp,ion,carboxy,brush,plasma,proj,adhesion,junction,lipo, ncol = 3, align = 'hv')

#paper figure 3a
save_plot("./figures/main/Fig3_all_clusters_blue.pdf", clusters_plot,base_height=15 )
```

Define further plot functions

``` r
plot_gene<-function(gene_selection){
ggplot(filter(zones,gene==gene_selection), aes(x=position, y=mean))+
  geom_rect(xmin = -Inf, xmax = 1.5,   ymin = -Inf, ymax = Inf,   fill = "#E2D9D9") +
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se), fill = brewer.pal(9,"Blues")[6]) +
  geom_line(color=brewer.pal(9,"Blues")[9],size=1.5)+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
  labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),axis.title.x=element_blank())+
  ylab("UMI fraction")+
  ggtitle(gene_selection)
}

plot_gene_small<-function(gene_selection){
  ggplot(filter(zones,gene==gene_selection), aes(x=position, y=mean))+
    geom_rect(xmin = -Inf, xmax = 1.5,   ymin = -Inf, ymax = Inf,   fill = brewer.pal(9,"Blues")[3]) +
    geom_ribbon(aes(ymin = mean-se, ymax = mean+se), fill = brewer.pal(9,"Blues")[6]) +
    geom_line(color=brewer.pal(9,"Blues")[9],size=1.5)+
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
    labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    ggtitle(gene_selection)
}
```

Plot zonation plots

``` r
#paper figure 6e
save_plot("./figures/main/Fig6e_enpp3.pdf", plot_gene_small("Enpp3"),base_height=2 )
save_plot("./figures/main/Fig6e_nt5e.pdf", plot_gene_small("Nt5e"),base_height=2 )
save_plot("./figures/main/Fig6e_ada.pdf", plot_gene_small("Ada"),base_height=2 )
save_plot("./figures/main/Fig6e_slc28a2.pdf", plot_gene_small("Slc28a2"),base_height=2 )

#paper figure 5
save_plot("./figures/main/Fig5b_Slc5a1.pdf", plot_gene("Slc5a1"),base_height=3 )
save_plot("./figures/main/Fig5d_Apob.pdf", plot_gene("Apob"),base_height=3 )

#paper figure 6a
save_plot("./figures/main/Fig6a_Egfr.pdf", plot_gene("Egfr"),base_height=3 )
save_plot("./figures/main/Fig6a_Klf4.pdf", plot_gene("Klf4"),base_height=3 )
save_plot("./figures/main/Fig6a_Fos.pdf", plot_gene("Fos"),base_height=3 )
save_plot("./figures/main/Fig6a_Junb.pdf", plot_gene("Junb"),base_height=3 )


#paper supp figures
save_plot("./figures/supplement/FigS5c_Nt5e.pdf", plot_gene("Nt5e"),base_height=3 )
save_plot("./figures/supplement/FigS4i_malat1.pdf", plot_gene("Malat1"),base_height=3 )
save_plot("./figures/supplement/FigS4G_neat1.pdf", plot_gene("Neat1"),base_height=3 )
save_plot("./figures/supplement/FigS5e_Slc28a2.pdf", plot_gene("Slc28a2"),base_height=3 )
save_plot("./figures/supplement/FigS5a_Cdh1.pdf", plot_gene("Cdh1"),base_height=3 )
save_plot("./figures/supplement/FigS4C_Gstm3.pdf", plot_gene("Gstm3"),base_height=3 )
save_plot("./figures/supplement/FigS4E_Nlrp6.pdf", plot_gene("Nlrp6"),base_height=3 )
save_plot("./figures/supplement/FigS2E_Reg1.pdf", plot_gene("Reg1"),base_height=3 )

save_plot("./figures/supplement/FigS7_Tfrc.pdf", plot_gene("Tfrc"),base_height=3 )
save_plot("./figures/supplement/FigS7_Reg3b.pdf", plot_gene("Reg3b"),base_height=3 )
save_plot("./figures/supplement/FigS7_Slc5a1.pdf", plot_gene("Slc5a1"),base_height=3 )
save_plot("./figures/supplement/FigS7_Cdh1.pdf", plot_gene("Cdh1"),base_height=3 )
save_plot("./figures/supplement/FigS7_Nt5e.pdf", plot_gene("Nt5e"),base_height=3 )
```

Scale zonation table to maximum and produce heatmaps

``` r
rownames(zonation)<-zonation$gene
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
A1 = zonation%>%
  select(2:8)
scaled<-A1/apply(A1,1,max)
scaled$gene<-rownames(scaled)
colnames(scaled) <- c(1,2,3,4,5,6,7,"gene")
scaled_clean<-na.omit(scaled)

scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%c("Reg3g","Reg3a","Reg3b","Lypd8","Nlrp6","Il18","Ccl25"))


cols<-c("#342D85","#2581C4","#35B4A2","#EAB94F","#EEE40B")

regs<-ggplot(scaled_zones, aes(position, gene)) + 
  geom_tile(aes(fill = mean), color = "black",size=0.5) +
  scale_fill_gradientn(colors=cols,guide = F)+
  coord_equal()+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line.x = element_blank(),
        axis.line.y = element_blank())

save_plot("./figures/main/Fig4a_regs.pdf", regs,base_width=6 )


scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%c("Slc2a2","Slc2a5","Slc5a1"))

carbo<-ggplot(scaled_zones, aes(position, gene)) + 
  geom_tile(aes(fill = mean), color = "black",size=0.5) +
  scale_fill_gradientn(colors=cols,guide = F)+
  coord_equal()+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line.x = element_blank(),
        axis.line.y = element_blank(),axis.text.x = element_blank())

scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%c("Slc15a1"))

protein<-ggplot(scaled_zones, aes(position, gene)) + 
  geom_tile(aes(fill = mean), color = "black",size=0.5) +
  scale_fill_gradientn(colors=cols,guide = F)+
  coord_equal()+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line.x = element_blank(),
        axis.line.y = element_blank(),axis.text.x = element_blank())

scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%c("Slc7a7","Slc7a8","Slc7a9"))

aa<-ggplot(scaled_zones, aes(position, gene)) + 
  geom_tile(aes(fill = mean), color = "black",size=0.5) +
  scale_fill_gradientn(colors=cols,guide = F)+
  coord_equal()+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line.x = element_blank(),
        axis.line.y = element_blank(),axis.text.x = element_blank())




scaled_zones<-gather(scaled,key=position,value=mean,-c(gene))%>%
  filter(gene%in%c("Apob","Apobec1","Apoa1","Apoa4","Npc1l1"))

apo<-ggplot(scaled_zones, aes(position, factor(gene,levels = c("Npc1l1","Apoa1","Apoa4","Apob","Apobec1")))) + 
  geom_tile(aes(fill = mean), color = "black",size=0.5) +
  scale_fill_gradientn(colors=cols,guide = F)+
  coord_equal()+
  scale_x_discrete(breaks=c(1,2,3,4,5,6,7),labels=c("Crypt","V1","V2", "V3","V4","V5", "V6"))+
  theme(legend.title=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line.x = element_blank(),
        axis.line.y = element_blank(),axis.text.x = element_blank())


save_plot("./figures/main/Fig5a_carbo.pdf", 
          carbo,base_width =3 )
save_plot("./figures/main/Fig5a_peptide.pdf", 
          protein,base_width =3 )
save_plot("./figures/main/Fig5a_apo.pdf", 
          apo,base_width =3 )
save_plot("./figures/main/Fig5a_aa.pdf", 
          aa,base_width =3 )
```

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bindrcpp_0.2.2     readr_1.1.1        tidyr_0.8.1       
    ## [4] seqinr_3.4-5       RColorBrewer_1.1-2 cowplot_0.9.3     
    ## [7] ggplot2_3.0.0      dplyr_0.7.6       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.18     pillar_1.3.0     compiler_3.5.0   plyr_1.8.4      
    ##  [5] bindr_0.1.1      tools_3.5.0      digest_0.6.15    packrat_0.4.9-3 
    ##  [9] evaluate_0.11    tibble_1.4.2     gtable_0.2.0     pkgconfig_2.0.1 
    ## [13] rlang_0.2.1      yaml_2.2.0       withr_2.1.2      stringr_1.3.1   
    ## [17] knitr_1.20       hms_0.4.2        rprojroot_1.3-2  ade4_1.7-11     
    ## [21] grid_3.5.0       tidyselect_0.2.4 glue_1.3.0       R6_2.2.2        
    ## [25] rmarkdown_1.10   purrr_0.2.5      magrittr_1.5     backports_1.1.2 
    ## [29] scales_0.5.0     htmltools_0.3.6  MASS_7.3-50      assertthat_0.2.0
    ## [33] colorspace_1.3-2 labeling_0.3     stringi_1.2.4    lazyeval_0.2.1  
    ## [37] munsell_0.5.0    crayon_1.3.4
