# H3K4 methylation chipseq webapp
 This is the webapp image hosted in pythonanywhere
 Unlike previous versions the data is no longer kept in memory but in a hdf5 file
 Data_load.py creates the hdf5 file from bed graph files.

This branch will contain major changes in the web app in order to make it easier to costumize,
as of January 2016 the only change that seems complete and stable is that Dataload.py now creates a user defined hdf5 file.
