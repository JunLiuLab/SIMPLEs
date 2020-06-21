
facs_objects = TM_facs_mat.rds TM_facs_metadata.csv
droplet_objects = TM_droplet_mat.rds TM_droplet_metadata.csv

cell_ontology:
	wget -O cell_ontology.obo https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo

facs_all:
	wget -O TM_facs_mat.rds https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_mat.rds
	wget -O TM_facs_metadata.csv https://raw.githubusercontent.com/czbiohub/tabula-muris-vignettes/master/data/TM_facs_metadata.csv

facs_meta_data:
	wget -O TM_facs_metadata.csv https://raw.githubusercontent.com/czbiohub/tabula-muris-vignettes/master/data/TM_facs_metadata.csv

facs_raw_data:
	wget -O TM_facs_mat.rds https://s3.amazonaws.com/czbiohub-tabula-muris/TM_facs_mat.rds

droplet_all:
	wget -O TM_droplet_mat.rds https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.rds
	wget -O TM_droplet_metadata.csv https://raw.githubusercontent.com/czbiohub/tabula-muris-vignettes/master/data/TM_droplet_metadata.csv

droplet_meta_data:
	wget -O TM_droplet_metadata.csv https://raw.githubusercontent.com/czbiohub/tabula-muris-vignettes/master/data/TM_droplet_metadata.csv

droplet_raw_data:
	wget -O TM_droplet_mat.rds https://s3.amazonaws.com/czbiohub-tabula-muris/TM_droplet_mat.rds


mv_facs:
	rm -rf ${facs_objects}

mv_droplet:
	rm -rf ${droplet_objects}
