#detach("package:RncMapping", unload=TRUE)
library(RncMapping)

connexion = cBioConnect()

# Selection of the study id
studies = listStudies(connexion, "leukemia")
st_id = studies$cancer_study_id[1]
ov_id = "ov_tcga_pub" # Published ovarian

#study_data = cBioStudy(st_id, genes_list="http://acsn.curie.fr/files/acsn_v1.0.gmt")
visualizer = cBioNCviz(ov_id, genes_list="file:///bioinfo/users/mdorel/ACSN_cBioportal_binding/test_data/genes_list", method="profiles")
#visualizer = cBioNCviz(ov_id, method="profiles")
#visualizer = cBioNCviz(ov_id, genes_list=c("MEK", "DNC"), method="genes")

saveInFiles(visualizer) # For NaviCell export
saveData(visualizer) # For easy sharing

