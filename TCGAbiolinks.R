library(TCGAbiolinks)

projects <- getGDCprojects()

query.met <- GDCquery(project=projects$project_id[4], legacy=TRUE, data.category="DNA methylation")
str(query.met)
str(projects)
projects$project_id[4]
projects$disease_type.0


View(query.met$results[[1]])
