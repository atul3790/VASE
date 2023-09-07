Keywords:

These are the different keywords you can use to get required observation tables from ALMA Query service. 
Eg: Keyword_ID=['science_keyword','band_list']   & Keyword_values=[['Main sequence stars','Brown dwarfs'],[3,6,9]]
For detailed list of values for different keywords visit: http://almascience.nrao.edu/aq/


Position
  Source name (Resolver)           : source_name_resolver               
  Source name (ALMA)               : source_name_alma                   
  RA Dec                           : ra_dec                             
  Galactic                         : galactic                           
  Angular resolution               : spatial_resolution                 
  Largest angular scale            : spatial_scale_max                  
  Field of view                    : fov                                

Energy
  Frequency                        : frequency                          
  Bandwidth                        : bandwidth                          
  Spectral resolution              : spectral_resolution                
  Band                             : band_list            -> 3(84-116 GHz)=3, 4(125-163 GHz)=4, 5(163-211 GHz)=5, 6(211-275 GHz)=6, 7(275-373 GHz)=7, 8(385-500 GHz)=8, 9(602-720 GHz)=9, 10(787-950 GHz)=10

Time
  Observation date                 : start_date                         
  Integration time                 : integration_time                   

Polarisation
  Polarisation type                : polarisation_type    -> Stokes I=0, Single=1, Dual=2, Full=3|4

Observation
  Line sensitivity (10 km/s)       : line_sensitivity                   
  Continuum sensitivity            : continuum_sensitivity              
  Water vapour                     : water_vapour                       

Project
  Project code                     : project_code                       
  Project title                    : project_title                      
  PI name                          : pi_name                            
  Proposal authors                 : proposal_authors                   
  Project abstract                 : project_abstract                   
  Publication count                : publication_count                  
  Science keyword                  : science_keyword      -> 

Publication
  Bibcode                          : bibcode                            
  Title                            : pub_title                          
  First author                     : first_author                       
  Authors                          : authors                            
  Abstract                         : pub_abstract                       
  Year                             : publication_year                   

Options
  (x) View:                        : result_view          = observation    
  ( ) View:                        : result_view          = project        
  ( ) View:                        : result_view          = biblio         
  [ ] public data only             : public_data          = public         
  [x] science observations only    : science_observations = =%TARGET%      


ALMA Query output table has the following columns. You may chose the ones relevant and supply it as a list 
Eg: ['Project code','Source name','RA',]

'Project code',
 'Source name',
 'RA',
 'Dec',
 'Galactic longitude',
 'Galactic latitude',
 'Band',
 'Spatial resolution',
 'Frequency resolution',
 'Array',
 'Mosaic',
 'Integration',
 'Release date',
 'Frequency support',
 'Velocity resolution',
 'Pol products',
 'Observation date',
 'PI name',
 'SB name',
 'Proposal authors',
 'Line sensitivity (10 km/s)',
 'Continuum sensitivity',
 'PWV',
 'Group ous id',
 'Member ous id',
 'Asdm uid',
 'Project title',
 'Project type',
 'Scan intent',
 'Field of view',
 'Largest angular scale',
 'QA2 Status',
 'COUNT',
 'Science keyword',
 'Scientific category',
 'ASA_PROJECT_CODE'
 
 ##################################################################################################################################
# Science Keywords usable to find the data for relevant science case (Following keywords only concern stellar / solar observations)
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Category 3 – ISM, star formation and astrochemistry

     Outflows, jets and ionized winds

     High-mass star formation

     Intermediate-mass star formation

     Low-mass star formation

     Pre-stellar cores, Infra-Red Dark Clouds (IRDC)

     Astrochemistry

     Inter-Stellar Medium (ISM)/Molecular clouds

     Photon-Dominated Regions (PDR)/X-Ray Dominated Regions (XDR)

     HII regions

     Magellanic Clouds

 
 Category 4 – Circumstellar disks, exoplanets and the solar system

     Debris disks

     Disks around low-mass stars

     Disks around high-mass stars

     Exoplanets

     Solar system: Comets

     Solar system: Planetary atmospheres

     Solar system: Planetary surfaces

     Solar system: Trans-Neptunian Objects (TNOs)

     Solar system: Asteroids


 Category 5 – Stellar evolution and the Sun

     The Sun

     Main sequence stars

     Asymptotic Giant Branch (AGB) stars

     Post-AGB stars

     Hypergiants

     Evolved stars: Shaping/physical structure

     Evolved stars: Chemistry

     Cataclysmic stars

     Luminous Blue Variables (LBV)

     White dwarfs

     Brown dwarfs

     Supernovae (SN) ejecta

     Pulsars and neutron stars

     Black holes

     Transients