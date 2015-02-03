detach("package:RncMapping", unload=TRUE)
library(RncMapping)

conn = cBioConnect()

# Selection of the study id
studies = listStudies(conn, "leukemia")
st_id = studies$cancer_study_id[1]
ov_id = "ov_tcga_pub" # Published ovarian

#study_data = cBioStudy(st_id, genes_list="http://acsn.curie.fr/files/acsn_v1.0.gmt")
visualizer = cBioNCviz(st_id, genes_list="file:///bioinfo/users/mdorel/ACSN_cBioportal_binding/test_data/genes_list", method="genes")
visualizer = cBioNCviz(ov_id, method="profiles")

saveInFiles(visualizer) # For NaviCell export
saveData(visualizer) # For easy sharing

