path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")

library(TCGAbiolinks)
library(dplyr)
library(DT)

projects <- getGDCprojects()
projects
query.meth <- GDCquery("TCGA-LUAD", legacy=FALSE, data.category="DNA Methylation")
clin <- GDCquery_clinic("TCGA-LUAD", type="clinical")
biospec <- GDCquery_clinic("TCGA-LUAD", type="Biospecimen")

GDCdownload(query.meth, directory=path.data)



View(query.meth$results[[1]])

clin$bcr_patient_barcode %in% substr(query.meth$results[[1]]$cases, 1, 12)


datatable(as.data.frame(colData(clin)), 
          options=list(scrollX=TRUE, keys=TRUE, pageLength=5), rownames = FALSE)

summary(query.clin)
summary(query.meth)

summary(query)
TCGAbiolinks:::getProjectSummary("TCGA-LUAD")

data <- GDCprepare(query, directory=path.data)
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)