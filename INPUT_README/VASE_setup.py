###############################
##### Directories & Input files
###############################

basedir='/mn/stornext/d19/RoCS/alma/atulm/Post_Doc_EMISSA/Codes/VASE'
Table_dir='/mn/stornext/d19/RoCS/alma/atulm/Post_Doc_EMISSA/Tables/VASE_tables'# Code will create this if not present
Prop_file= 'Vizier_prop_keywords.txt'# What properties you need? Fill up the file Vizier_prop_keywords.txt saved in basedir. B-V color is must
Catalog_pref_file='Imp_catalogs_Vizierlink.txt' #Preferred catalogues to be looked up first. Add Catalogs to Imp_catalogs_Vizierlink.txt in basedir
stars_file='Stars_list.txt' #Stars list
INPUT_files_folder='INPUT_README' # Must place this within the basedir. The location where input txt files are present.
email='atul.3.7.1990@gmail.com' #EMAIL to send the work finished notice.
Option=3 #'1. Make star catalog 2. Download Catalogs 3. Search for and catalog ALMA observations of interest using keywords 4.Find & download all ALMA data available for set of objects. 

#### GENERAL QUERY OUTPUT SETTING ###########
Row_limit=999999999 #Limit on numer of rows in a downloaded catalog. Set this really large if the purpose is to download entire catalog
Waittime=240 # Time in seconds to wait till the website responds to your query. Poor connectivity and stars with lot of observations would require more time for Vizier server to collect all data and respond to user query. 

Buffer_cpus=10 # Total number of cores to be left behind while parallel processing is done.
Nice=5 # CPU nice priority value. lower the value more priority your processes get. default unix number is 0. range : -20 to 19
###############################
##### INPUT for Option 1
###############################

look_across_cats=True #Do you need to include other catalogs in the search for info apart from Preferred ones.?
radius_max=3 # Search for object within this max radius in arssec
radius_min=1.# Minimum radius of the region within which the source should be looked up for.
MY_catalog='Test_opt1'# Name of the new catalog to be made. Give just the name and no extensions as .txt, .csv etc. 
							#Final output will be .tsv readable in spreadsheets.
Only_make_catalog=False #If you have already saved star data as pickle files.
is_new_catalog=True #If you need to make a whole new catalog tsv file. Else code will look for existing tsv file and append the new stars to it.
delete_old=False # If you want to delete any catalog with same name existing the catalog directory
append_existing_logs=False # Append logs in existing logfiles
delete_old_logs=False
Calc_params=True # Do you want the code to calculate values of unknown parameters using known stellar correlation curves.
				##WORKS ONLY FOR MAIN SEQUENCE STARS
val_prec=4 #Decimal places to which parameters have to be calculated and reported. Mind this only if Calc_params=True				
update_bib=False #If reference txt file is to be updated with new references.
bib_txtfile='Test_opt1_ref.txt' #If update_bib=True then give the name of the reference file.

###############################
##### INPUT for Option 2
###############################

Catalog_names_file='Catalogs_ToDownload.txt' # Text file containing names of catalogs to be downloaded and saved.
Catalog_dir='Useful_catalogs' # Location to save catalogs will be made inside Table_dir
 
###############################
##### INPUT for Option 3
###############################

#For all keywords and relevant IDs, and for ALMA table column names please look at 'ALMA_Query_help_README.txt'
Keyword_IDs_list=[['science_keyword']]*2 # Search keywords that helps to find relevant ALMA data. This is a must entry. List of lists. Each list will be independently used to seacrh for observations matching the criteria specified in the list
Keyword_values_list=[['Luminous Blue Variables (LBV)'],['Main sequence stars']] # Value for the keywords you choose above. 
										#Eg: Keyword_IDs_list=[['science_keyword'],['source_name_alma']] ; Keyword_values_list=[['Main sequence stars'],['The Sun']]. 2 catalogs one for Main sequence stars and 1 for Sun will be made.
catalog_names=['LBV_Aug2022','MS_Aug2022'] # Just give the name without any extension as a list. Eg [xxx] if only 1 keyword is input.
ALMA_catalog_dir='ALMA_Catalogs' # This will be made inside Table_dir if the folder isn't present already.
columns_wanted=['Project code','Observation date','Source name','RA','Dec', 'Galactic longitude','Galactic latitude', 'Field of view','Band',
'Spatial resolution','Frequency resolution','Integration','Frequency support', 'Pol products'] # List of column names you need to have in the final saved catalog. LEave blank if you need all columns
SIMBAD_ID_needed=True
imp_props=['Band','Integration','Frequency support','Field of view']# Give the column names / Kerword_IDs for important properties which needs to be appended to the row data corresponding 
			#to repeated stars. One row will be made for each star in te table by appending the values of these select properties for multiple
			#observations of the same star. Project code is anyway default element in imp_props.					
Public=True #Return only publicly available ALMA obs Metadata
Science_only=True # Return only data marked science
##### Repeater objects filtering######
filter_repeaters=True # Some stars are repeatedly observed and recorded in ALMA catalog. If this is set true only the first observation of the object will be displayed in table. The code will record the project code of any later observation and put it in the final table. 
# If filter_repeaters == True. The following variable value is relevant. Else neglect it
known_alias_names=[['Epsilon_Eridani','Eps_Eri','Epsilon Eri'],['Alpha Centaury','Alpha cen'],['Beta_Pictoris','beta pic']]# Works if filter_repeater=True
filter_repeats_with_SIMBAD_info=True #Use this to filter to eliminate multiple entries of an object with multiple names.
							#!!!!!! Dont use this for cataloging individual sources like Sun. 

#################################
##### INPUT for Option 4
#################################
url='https://jvo.nao.ac.jp/portal/alma/alma-meta.psv.gz' # File with meta data of all ALMA observations so far
All_almaobs_file='ALL_ALMA_obs_sofar' # Just give the name of the file.File is expected to be present in ALMA_catalog_dir
Search_stars_meta=False # If you just need to download the ALMA metadata file for search stars and tabulate.
Download_all=False #Do you need the code auto download all fits data. This works if Search_stars_meta==False
just_get_SearchStarimgs=True # If metadata is already downloaded and now just FITS files are to be downloaded.
Get_only_Quicklookimg=True # If True fits images wont be downloaded and only JVO spectra and jpg images will be saved. JVO images however will be saved if this is set False.
Search_stars='ALMA_Search_list.txt' # Stars or objects to be search for in ALMA observation catalog.
Search_result_dir='ALMA_search_results' # will be saved in Table_dir/MY_tables
Out_table_name='A-Table_MS_stars_B36'# Output details of search stars
my_dpi=96 # jpg image dpi.
Is_SIMBAD_names=True # Are the names provided SIMBAD IDs?
ALMA_obs_cat_pfile='ALL_ALMA_stellar_obssofar' #Table pickle file to look up for ALMA observations for the given list of objects. If this pickle file isnt generated, the code will generate it.
rel_otype=['*'] # Relevant SIMBAD Otype for the objects you need the code to record in the pfile. If these are stars use '*'. No need to give specific type names. Leave it as [] if you need the code to catalog details of all objects in ALMA data repo.
strict_filter=False # What if simbad query failed to give an otype? should you still include the object in the file (strict_filet=True) or ignore it?
