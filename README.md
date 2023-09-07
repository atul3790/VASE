# VASE
VizieR ALMA and SIMBAD Explorer
##################### Vizier ALMA and SIMBAD Explorer (VASE) Manual ###################### ##########################################################################################

VASE is a package which can be used to retrieve information about any astronomical object. It is however specially tuned for stellar observations. The directory structure of VASE is shown below.

Folders:

INPUT_README
------------
1. VASE_setup.py : parameter file with detailed choice list. Update this as per the needs before running VASE
2. Stars_list.txt: List of stars to search for in VizieR and SIMBAD database
3. ALMA_Search_list.txt: List of stars to search for in ALMA database

4. Catalogs_ToDownload.txt: Names of catalogs to download with the ViZieR folder path (eg: B/230)
5. Imp_catalogs_Vizierlink.txt: Names and folder paths of catalogs which are important while searching for sources
6. Property_keywords_allowed.txt: Keywords of stellar properties which can be querried - USe as a initial successfully tried list. Update things with time.
7. Vizier_prop_keywords.txt: 3 column text file with property names, list of variable names used in various catalogs to refer to the quantity, data type and UCD - a code VizieR gives to variables to specify its physical meaning.
8. ALMA_Query_help_README.txt: different keywords to use while searching for ALMA data.

Logs
----
After each run, important information about succesful and failed searches will be reported in this directory.
