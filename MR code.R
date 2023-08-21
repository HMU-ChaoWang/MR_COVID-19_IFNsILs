library(remotes)
library(TwoSampleMR)
#IL
ILall <- read.csv("Interleukin.csv",sep = ",")
ILnum <- ILall$GWAS.ID
ILname <- ILall$Trait

for (i in 1:length(ILname)) {
  dir.create(paste0("E:/MR/IL/1e-5/",ILname[i]))
}

for (i in 1:length(ILnum)) {
  IL <- extract_instruments(
    outcomes = ILnum[i],
    clump = T,p1=1e-5,
    r2=0.01,kb=1,
    access_token = NULL
  )
  
  COVID19 <- extract_outcome_data(
    snps = IL$SNP,
    outcomes = "ebi-a-GCST011074",
    proxies = F,
    maf_threshold = 0.01,
    access_token = NULL
  )
  data <- harmonise_data(exposure_dat = IL,outcome_dat = COVID19,action = 2)
  write.csv(data,paste0("E:/MR/IL/1e-5/",ILname[i],"_data.csv"))
  res <- mr(data)
  write.csv(res,paste0("E:/MR/IL/1e-5/",ILname[i],"_res.csv"))
  #异质性检验
  heterogeneity <- mr_heterogeneity(data)
  write.csv(heterogeneity,paste0("E:/MR/IL/1e-5/",ILname[i],"_heterogeneity.csv"))

  pleiotropy <- mr_pleiotropy_test(data)
  write.csv(pleiotropy,paste0("E:/MR/IL/1e-5/",ILname[i],"_pleiotropy.csv"))

  single1<-mr_leaveoneout(data)
  res_single1<-mr_singlesnp(data)
  
  pdf(paste0("E:/MR/IL/1e-5/",ILname[i],"_leaveplot.pdf"))
  leaveplot<- mr_leaveoneout_plot(single1) 
  print(leaveplot)
  dev.off()
  
  pdf(paste0("E:/MR/IL/1e-5/",ILname[i],"_scatter.pdf"))
  scatter <- mr_scatter_plot(res,data)
  print(scatter)
  dev.off()
  
  pdf(paste0("E:/MR/IL/1e-5/",ILname[i],"_forestplot.pdf"))
  forestplot <- mr_forest_plot(res_single1)
  print(forestplot)
  dev.off()
  
  #ggsave(paste0("E:/MR/IL/1e-5/",ILname[i],"_allplot.pdf"), plot = allplot, width = 12, height = 10)
  ORratio <- generate_odds_ratios(res)
  write.csv(ORratio,paste0("E:/MR/IL/1e-5/",ILname[i],"_ORratio.csv"))
  print(i)
}
  
#IFN
setwd("E:/MR/IFN/")
IFNall <- read.csv("IFN.csv",sep = ",")
IFNnum <- IFNall$GWAS.ID
IFNname <- IFNall$Trait
length(IFNnum)

length(IFNnum)
i=10

for (i in 1:length(IFNnum)) {
  IFN <- extract_instruments(
    outcomes = IFNnum[i],
    clump = T,p1=1e-5,
    r2=0.01,kb=1,
    access_token = NULL
  )
  
  COVID19 <- extract_outcome_data(
    snps = IFN$SNP,
    outcomes = "ebi-a-GCST011074",
    proxies = F,
    maf_threshold = 0.01,
    access_token = NULL
  )
  
  data <- harmonise_data(exposure_dat = IFN,outcome_dat = COVID19,action = 2)
  write.csv(data,paste0("E:/MR/IFN/1e-5/",IFNname[i],"_data.csv"))
  res <- mr(data)
  write.csv(res,paste0("E:/MR/IFN/1e-5/",IFNname[i],"_res.csv"))
  
  heterogeneity <- mr_heterogeneity(data)
  write.csv(heterogeneity,paste0("E:/MR/IFN/1e-5/",IFNname[i],"_heterogeneity.csv"))

  pleiotropy <- mr_pleiotropy_test(data)
  write.csv(pleiotropy,paste0("E:/MR/IFN/1e-5/",IFNname[i],"_pleiotropy.csv"))

  single1<-mr_leaveoneout(data)
  leaveplot<- mr_leaveoneout_plot(single1) 
  pdf(paste0("E:/MR/IFN/1e-5/",IFNname[i],"_leaveplot.pdf"))
  leaveplot
  dev.off()
  scatter <- mr_scatter_plot(res,data)
  pdf(paste0("E:/MR/IFN/1e-5/",IFNname[i],"_scatter.pdf"))
  scatter
  dev.off()

  forestplot <- mr_forest_plot(res_single1)
  pdf(paste0("E:/MR/IFN/1e-5/",IFNname[i],"_forestplot.pdf"))
  forestplot
  dev.off()
  #ggsave(paste0("E:/MR/IFN/1e-5/",IFNname[i],"_allplot.pdf"), plot = allplot, width = 12, height = 10)
  ORratio <- generate_odds_ratios(res)
  write.csv(ORratio,paste0("E:/MR/IFN/1e-5/",IFNname[i],"_ORratio.csv"))
  print(i)
}
