
conn = CGDS("http://www.cbioportal.org/")

# Selection of the study id
studies = getCancerStudies(conn)
st_id = studies$cancer_study_id[1]

# Select all patients
#cases = getCaseLists(conn, st_id) # Use and display information ?
ca_id = paste0(st_id, "_all")
# Retrieve genetic profiles ids of all profiles
profiles = getGeneticProfiles(conn, st_id)
pr_id = profiles$genetic_profile_id

genes_data = retrieveDataSet(conn, pr_id, ca_id)

