library(RncMapping)

conn = cBioConnect()

# Selection of the study id
studies = listStudies(conn)
st_id = studies$cancer_study_id[1]

genes_data = cBioStudy(st_id)

visualizer = NCviz(nc_url="http://acsn.curie.fr/files/acsn_v1.0.owl", cell_type="Acute Myeloid Leukemia", cbio_data=genes_data)

