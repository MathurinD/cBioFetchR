library(RncMapping)

conn = cBioConnect()

# Selection of the study id
studies = listStudies(conn)
st_id = studies$cancer_study_id[1]

#visualizer = cBioStudy(st_id, genes_url="http://acsn.curie.fr/files/acsn_v1.0.gmt")
visualizer = cBioNCviz(st_id, genes_url="file:///bioinfo/users/mdorel/ACSN_cBioportal_binding/test_data/genes_list", method="genes")

saveInFiles(visualizer)

