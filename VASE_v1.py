'''
Code gathers information about various properties of stars by browsing through catalogues in Vizier database.
User need to give 3 text files as Input.

'''
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.table import Table,Column
import pickle,importlib,urllib,requests
from multiprocessing import Pool
import numpy as np
import astropy.units as u
import os,ast,sys
import time,glob
from astropy.io import fits
import matplotlib.pyplot as plt 
from pandas import unique
from astroquery.alma import Alma
from Star_spec import *

########### INPUT ######################################################################
import VASE_setup
importlib.reload(VASE_setup)
from VASE_setup import *
#########################################################################################

#Session=requests.Session()

## Checking int or float for numbers##########################
def is_float(num):
	tp=0
	try:

		tp=float(num)
		if ~np.isnan(tp):
			return True
		else:
			return False	
	except:
		return False
					
def is_int(num):
	if str(num).count('.')==0:
		try:
			tp=int(num)
			if ~np.isnan(tp):
				return True
			else:
				return False	
		except:
			return False
	else:
		return False


##################### INITIAL SETUP of QUERY modules ##############	
os.chdir(basedir)
Vizier.ROW_LIMIT=Row_limit
Vizier.TIMEOUT=Waittime
Alma.TIMEOUT=Waittime
Simbad.TIMEOUT=Waittime
Simbad.add_votable_fields('ra','dec','otypes','sp','otype','flux(B)','flux(V)','distance','diameter')

################## SYSTEM Details #################################################
compd=os.uname().nodename
Cores=os.cpu_count()
Nproc=np.min([Cores-2,Cores-Buffer_cpus])
jvo_imgurl="https://jvo.nao.ac.jp/portal/alma/archive.do?pictSize=512&dataType=image&action=quicklook&dataId="
jvo_specurl="https://jvo.nao.ac.jp/portal/alma/archive.do?pictSize=512&dataType=spect&action=quicklook&dataId="
############# Polishing Input files to be read by functions ##########################
ckeys,cvals=np.loadtxt(INPUT_files_folder+'/'+Catalog_pref_file,dtype=str,unpack=1,usecols=(0,1),delimiter='\t')
Catalog_pref=dict(zip(ckeys,cvals))
ckeys,cvals,typs,ucds=np.loadtxt(INPUT_files_folder+'/'+Prop_file,dtype=str,unpack=1,usecols=(0,1,2,3),delimiter='\t\t')
cvals=[ast.literal_eval(i) for i in cvals]
ucds=[ast.literal_eval(i) for i in ucds]
tcvals=[]
tucds=[]
key_typ=dict(zip(ckeys,typs))
for ele in range(len(cvals)):
	tcvals+=cvals[ele]
	tucds+=ucds[ele]
tcvals=[i.upper() for i in tcvals]	
tucds1=','.join(tucds).upper()
#print(ckeys,'\n',tcvals,'\n',tucds)	
if 'B-V' in ckeys or 'B-V' in tcvals or 'B-V' in tucds1:
	pass
else:	
	ckeys+=['B-V']
	cvals+=[['B-V']]
	typs+=['float']
	ucds+=[['PHOT_CI_B-V','PHOT_JHN_B-V']]
del tucds
del tcvals
del tucds1
if Calc_params==True:
	if 'Distance' not in ckeys:
		ckeys+=['Distance']
		cvals+=[['Dist']]
		typs+=['float']
		ucds+=[['POS_GAL_HC']]
	if 'M' not in ckeys:
		ckeys+=['M']
		cvals+=[['Mass','Mstar','M']]
		typs+=['float']
		ucds+=[['PHYS_MASS_MISC']]
	if 'Teff' not in ckeys:
		ckeys+=['Teff']
		cvals+=[['Teff']]
		typs+=['float']
		ucds+=[['PHYS_TEMP_EFFEC','PHYS_TEMP_MISC']]
	if 'Lbol' not in ckeys:
		ckeys+=['Lbol']
		cvals+=[['L','LBol','Lum','Lstar']]
		typs+=['float']
		ucds+=[['PHOT_BOL_LUMINOSITY','PHYS_LUMINOSITY_GENERAL']]
	if 'logRHK' not in ckeys:
		ckeys+=['logRHK']
		cvals+=[['logRHK','log10RHK']]
		typs+=['float']
		ucds+=[['SPECT_INDEX_MISC','PHYS_EM-MEASURE','PHOT_FLUX_RATIO','SPECT_FLUX_RATIO']]
	if 'S' not in ckeys:
		ckeys+=['S']
		cvals+=[['S','GSval','SMW']]
		typs+=['float']
		ucds+=[['SPECT_INDEX_MISC','PHOT_FLUX_RATIO','PHYS_EM-MEASURE','SPECT_FLUX_RATIO']]	
# Assigning data type checking functions to each property  
typs1=[is_float if i == 'float' else str.isalpha if i=='str' else is_int if i=='int' else str.isalnum for i in typs ]
cvals1=[[cvals[i],typs1[i],ucds[i]] for i in range(len(cvals))]
Props=dict(zip(ckeys,cvals1))

stars=np.loadtxt(INPUT_files_folder+'/'+stars_file,dtype=str,delimiter=';')

imp_props=list(np.unique(np.array(['Observation date']+imp_props))) 
if 'Frequency support' not in imp_props:
	imp_props+=['Frequency support']
if 'Frequency resolution' not in imp_props:
	imp_props+=['Frequency resolution']
if 'Pol products' not in imp_props:			
	imp_props+=['Pol products']
imp_props+=['Frequency range']

for i in range(len(known_alias_names)):
	tmp=[j.replace(' ','').replace('_','').upper() for j in known_alias_names[i]]
	known_alias_names[i]=tmp
columns_wanted = ['Observation date']+columns_wanted if 'Observation date' not in columns_wanted else columns_wanted
	
########### Checking directory structure ###################
if not os.path.isdir(Table_dir):
	os.mkdir(Table_dir)
Catalog_dir=Table_dir+'/'+Catalog_dir+'/'
if Option==2:
	if not os.path.isdir(Catalog_dir):
		os.mkdir(Catalog_dir)
ALMA_catalog_dir=Table_dir+'/'+ALMA_catalog_dir
if Option==3:
	if not os.path.isdir(ALMA_catalog_dir):
		os.mkdir(ALMA_catalog_dir)
if Option==1:
	Table_dir+='/MY_tables/'
	if not os.path.isdir(Table_dir):
		os.mkdir(Table_dir)							
	if not os.path.isdir(Table_dir+'Data/'):
		os.mkdir(Table_dir+'Data/')							
	pfildir=Table_dir+'Data/'+MY_catalog+'_Stardata/'
	if not os.path.isdir(pfildir):
		os.system('mkdir '+pfildir)
	pfildir_cats=pfildir+MY_catalog+'_Catalogs/'
	if not os.path.isdir(pfildir_cats):
		os.system('mkdir '+pfildir_cats)
'''		
#####################################################################################################		
######## Dealing with ALMA obs search separately ###################################################
#####################################################################################################	
'''
if Option==4:
	plt.ioff()
	###################### Initialisation of variables ##################		
	search_stars=np.loadtxt(INPUT_files_folder+'/'+Search_stars,dtype=str,delimiter=';')
	FLAG=0		
	Srcs=[]
	goods=[]
	src_list=np.array([])
	src_uniq_list=np.array([])
	Goodps=[]
	Good_srcs=[]
	Good_IDs=[]
	Badps=[]
	Bad_noIDs=[]
	Bad_withIDs=[]
	Bad_IDs=[]
	Allps=[]
	Good_dets=[]
	Bad_dets=[]
	Bad_noIDdets=[]
	allID_match=dict(zip([],[]))
	allGoodpfs=[]
	allGoodIDs=np.array([])
	allBadpfs=[]
	allBadIDs=np.array([])
	id_matcher=dict(zip([],[]))
	src_matcher=dict(zip([],[]))
	res_fold_matcher=dict(zip([],[]))
#############################################################################
###### Folder stucture making/ensuring ######################################
	
	if not os.path.isdir(Table_dir):
		os.mkdir(Table_dir)
	if not os.path.isdir(ALMA_catalog_dir):
		print('No Catalogs available in the user specified folder, ',ALMA_catalog_dir)
	os.chdir(ALMA_catalog_dir)
	SR_dir=Table_dir+'/MY_tables/'+Search_result_dir+'/'
	if not os.path.isdir(SR_dir):
		os.mkdir(SR_dir)
	Outtabdir=SR_dir+Out_table_name+'/'
	if not os.path.isdir(Outtabdir):
		os.system('mkdir '+Outtabdir)
	if not os.path.isfile(All_almaobs_file+'.psv'):
		urllib.request.urlretrieve('https://jvo.nao.ac.jp/portal/alma/archive.do?action=download.metadata',All_almaobs_file+'.psv.gz')
		os.system('gunzip '+All_almaobs_file+'.psv.gz')
	if not os.path.isdir(All_almaobs_file+'_pfils/'):
		FLAG=1
	elif len(glob.glob(All_almaobs_file+'_pfils/*.p'))==0:
		FLAG=1
	if FLAG==1:	
		allf=np.genfromtxt(All_almaobs_file+'.psv',delimiter='\n',dtype=str)
		lens=np.array([len(i.split('|')) for i in allf])
		lens_uniq=unique(lens)
		maxl=0
		max_val=0
		for i in lens_uniq:
			if len(np.where(lens==i)[0])>maxl:
				max_val=np.where(lens==i)[0]
				maxl=len(max_val)
		goods=allf[max_val]
		src_list=np.array([i.split('|')[4] for i in goods])
		src_uniq_list=unique(src_list)
		#goods=[i.split('|') for i in goods]
			
		if os.path.isdir(All_almaobs_file+'_pfils/'): # If the complete ALMA observation record p file isnt there then better make from scratch
			os.system('rm -rf '+All_almaobs_file+'_pfils/')
		os.mkdir(All_almaobs_file+'_pfils/')
		Srcs=[src_uniq_list[i]+';;Mastercatalog;;'+str(i) for i in range(len(src_uniq_list))]
		ALMA_pfilmatch=dict(zip(['Mastercatalog','Subcatalog'],[All_almaobs_file+'_pfils/',ALMA_obs_cat_pfile+'_pfils/']))		

		fil=open('ALL_ALMA_obs.txt','a')
		fil.write('#ALMA name\tSIMBAD ID\tOTypes\tLoc in Sourcelist')
		fil.close()	
		
	########################################################################################		
	'''
	######################################################################################			
	###### Finding SIMBAD IDS for sources ################################################
	######################################################################################
	'''

	def get_SIMBAD_ids(idl):
		idn=idl.split(';;')[0]
		pfildir=ALMA_pfilmatch[idl.split(';;')[1]]
		locinsl=idl.split(';;')[2]
		filtag=idn.replace('/','_slh_')
		Rel_Otypes=[]
		Strict_filter=False
		if not pfildir==All_almaobs_file+'_pfils/':
			Rel_Otypes+=rel_otypes
			Strict_filter=False or strict_filter
			print('Strict_filter: ',Strict_filter)
		simID='--'

		fil=open('ALL_ALMA_obs.txt','a')	
		try:
			simoj=Simbad.query_object(idn)
			if 'NONE' in str(type(simoj)).upper() or 'UNRECOGNIZED' in str(type(simoj)).upper():
				print('No output for object, ', idn)
				filtag+='$_$noSIMRESULT'
				if Strict_filter==False:
					fil.write('\n'+idn+'\t--\t--\t'+locinsl)
					fil.close()
					pickle.dump((idn,'--','--',locinsl),open(pfildir+filtag+'$_$'+locinsl+'.p','wb'))
					return 0
				else:	
					print('Strictly filtered')
					return 0
			otypes=simoj['OTYPES'].data.data[0]
			simID=simoj['MAIN_ID'].data.data[0]
			if len(otypes)==0 and len(simID)==0:
				print('No output for object, ', idn)
				filtag+='$_$noSIMRESULT'
				if Strict_filter==False:
					fil.write('\n'+idn+'\t--\t--\t'+locinsl)
					fil.close()
					pickle.dump((idn,'--','--',locinsl),open(pfildir+filtag+'$_$'+locinsl+'.p','wb'))
					return 0
				else:	
					print('Strictly filtered')
					return 0
				
			if 'BYTES' in str(type(otypes)).upper() and len(otypes)>0:
				otypes=otypes.decode('UTF-8')
				otypes=otypes.split('|')
			elif 'STR' in str(type(otypes)).upper() and len(otypes)>0:
				otypes=otypes.split('|')
			else:
				otypes=['--']
				print('Bad OTYPE found!!!')
				filtag+='$_$noOTYPES'

			if len(Rel_Otypes)>0:
				if otypes==['--']:
					Rel_Otypes+=['--']
					if Strict_filter==True:
						print('Being strict, not recording the data')
						return 0
				if not any([i in otypes for i in Rel_Otypes]):
					print('No OType match the user requirements..')
					return 0				


			if 'BYTES' in str(type(simID)).upper() and len(simID)>0:
				simID=simID.decode('UTF-8')
			elif 'STR' in str(type(simID)).upper() and len(simID)>0:
				pass
			else:
				print('Bad ID type found!!!')
				simID='--'
				filtag+='$_$noSIMBADID'
			fil.write('\n'+idn+'\t'+simID+'\t'+';'.join(otypes)+'\t'+locinsl)
			fil.close()				
			pickle.dump((idn,simID,';'.join(otypes)),open(pfildir+filtag+'$_$'+locinsl+'.p','wb'))
		except:
			filtag+='$_$noQueryRESULT'
			print('Bad source ',idl)
			if Strict_filter==False:
				fil.write('\n'+idn+'\t--\t--\t'+locinsl)
				fil.close()
				pickle.dump((idn,'--','--'),open(pfildir+filtag+'$_$'+locinsl+'.p','wb'))	
		return 0		
	###############################################################################################################	

	if FLAG==1:
		print('Collecting SIMBAD IDs of all ALMA observed sources ad making pfils..')
		p=Pool(Nproc)
		os.nice(Nice)
		succ=p.map(get_SIMBAD_ids,Srcs)
		os.chdir(All_almaobs_file+'_pfils/')
		os.system('mkdir Good_ones/')
		os.system('mkdir Bad_ones_withSIMBADID/')
		os.system('mkdir Bad_ones_withnoSIMBADID/')
		
		########################## Analysing pickle files and segregating ones with and without SIMBAD ID/type ####

		Allps=glob.glob('*.p')
			
		for i in Allps:
			temV=pickle.load(open(i,'rb'))
			if not '$_$no' in i:
				Goodps+=[temV]
				Good_srcs+=[temV[0]]
				Good_IDs+=[temV[1]]
			else:
				Badps+=[temV]
				if not temV[1]=='--':
					Bad_withIDs+=[temV[0]]
					Bad_IDs+=[temV[1]]
				else:
					Bad_noIDs+=[temV[0]]
		Good_IDs=np.array(Good_IDs)
		Good_srcs=np.array(Good_srcs)
		Bad_IDs=np.array(Bad_IDs)
		Bad_noIDs=np.array(Bad_noIDs)
		Bad_withIDs=np.array(Bad_withIDs)		
		Bad_IDs_uniq=unique(Bad_IDs)
		Good_IDs_uniq=unique(Good_IDs)
		id_matcher=dict(zip(['good','bad'],[Good_IDs,Bad_IDs]))
		src_matcher=dict(zip(['good','bad'],[Good_srcs,Bad_withIDs]))
		res_fold_matcher=dict(zip(['good','bad'],['Good_ones/','Bad_ones_withSIMBADID/']))
				
		###########################################################################################################				
	'''
	#####################################################################################################
	###### Make pickle files for ALMA obs sources with SIMBAD ID #####################
	#####################################################################################################			
	'''	
	def pickle_sources_withID(src):
		i,tag=src.split(';;')
		Src_ar=src_matcher[tag]
		ID_ar=id_matcher[tag]
		alma_srcs=Src_ar[np.where(ID_ar==i)[0]]
		t_srcn=','.join(alma_srcs)
		t_ALMAID=[]
		t_obsdt=[]
		t_bnd=[]
		t_science=[]
		t_url=[]
		t_iscal=[]
		t_filname=[]
		locs=[]
		print('Analysing ',i)
		for Name in alma_srcs:
			locs+=list(np.where(src_list==Name)[0])
		for j in goods[locs]:
			tvls=j.split('|')
			t_ALMAID+=[tvls[0]]
			t_obsdt+=[tvls[6]]
			t_bnd+=[tvls[9]]
			t_science+=[tvls[2]]
			t_url+=[tvls[-1]]
			t_iscal+=[tvls[-2]]
			t_filname+=[tvls[-5]]
		Dets=[i,t_srcn,t_ALMAID,t_obsdt,t_bnd,t_science,t_url,t_iscal,t_filname] # Recording SIMBAD ID, ALMA Source name,ALMA dataID, Obs datetime, obs band, science case,fits file location.				
		pickle.dump(Dets,open(res_fold_matcher[tag]+i.replace('/','_slh_')+'.p','wb'))
		return 0	
	'''
	#####################################################################################################
	###### Make pickle files for ALMA obs sources with NO SIMBAD ID #####################
	#####################################################################################################			
	'''	

	def pickle_sources_noID(src):
		i=src.split('$_$')[0]
		locs=list(np.where(src_list==i)[0])
		t_srcn=i
		t_ALMAID=[]
		t_obsdt=[]
		t_bnd=[]
		t_science=[]
		t_url=[]
		locs=[]
		t_filname=[]
		t_iscal=[]
		print('Analysing ',i)
		for j in goods[locs]:
			tvls=j.split('|')
			t_ALMAID+=[tvls[0]]
			t_obsdt+=[tvls[6]]
			t_bnd+=[tvls[9]]
			t_science+=[tvls[2]]
			t_url+=[tvls[-1]]
			t_iscal+=[tvls[-2]]
			t_filname+=[tvls[-5]]
		Bad_noIDdets=[i,t_srcn,t_ALMAID,t_obsdt,t_bnd,t_science,t_url,t_iscal,t_filname] # Recording SIMBAD ID, ALMA Source name,ALMA dataID, Obs datetime, obs band, science case,fits file location.
		pickle.dump(Bad_noIDdets,open('Bad_ones_withnoSIMBADID/'+src.replace('/','_slh_')+'.p','wb'))				
		return 0				
	
	########################################################################################		
	if FLAG ==1:
		p=Pool(Nproc)
		os.nice(Nice)
		Gu=[i+';;good' for i in Good_IDs_uniq]
		Bu=[i+';;bad' for i in Bad_IDs_uniq]
		succ=p.map(pickle_sources_withID,Gu)
		print('Done with sources with SIMBAD entry...\n Moving on to the ones with SIMBAD ID but no Otype defined!!!')	

		succ=p.map(pickle_sources_withID,Bu)
		print('Done with ones with Otype undefined but SIMBAD ID defined. \n Moving on to ones with failed SIMBAD ID.')	

		Bad_noIDs_t=[Bad_noIDs[i]+'$_$'+str(i) for i in range(len(Bad_noIDs))]

		succ=p.map(pickle_sources_noID,Bad_noIDs_t)
		print('Done with recording all ALMA observed sources with no SIMBAD entry.')
		print('Good srcs: ',len(Good_IDs_uniq))
		print('Bad srcs with SIMBAD: ',len(Bad_IDs_uniq))
		print('Bad sources with no SIMBAD: ',len(Bad_noIDs))
		
	if Search_stars_meta==True:
		os.chdir(ALMA_catalog_dir+'/'+All_almaobs_file+'_pfils/')
		allGoodpfs=glob.glob('Good_ones/*.p')
		allGoodIDs=np.array([allGoodpfs[i].split('/')[1].split('$_$')[0].replace('_slh_','/').replace('.p','').replace(' ','')+'@_'+str(i) for i in range(len(allGoodpfs))])
		allBadpfs=glob.glob('Bad_ones_withSIMBADID/*.p')
		allBadIDs=np.array([allBadpfs[i].split('/')[1].split('$_$')[0].replace('_slh_','/').replace('.p','').replace(' ','')+'@_'+str(i) for i in range(len(allBadpfs))])
		Bad_noIDpfs=glob.glob('Bad_ones_withnoSIMBADID/*.p')
		Bad_noIDs_nosp=[i.split('/')[-1].split('$_$')[0].replace(' ','') for i in Bad_noIDs]
		allworstpfmatch=dict(zip(Bad_noIDs_nosp,Bad_noIDpfs))
		allID_match=dict(zip(['Good_ones/','Bad_ones_withSIMBADID/'],[allGoodIDs,allBadIDs]))
		allpf_match=dict(zip(['Good_ones/','Bad_ones_withSIMBADID/'],[allGoodpfs,allBadpfs]))
		print('No. of allGoodIDs: ',len(allGoodIDs))
		print('No. of allBadIDs: ',len(allBadIDs))
		
#########################################################################################		
	'''
	#############################################################################################################	
	#### Download ALMA metadata of the sources required and tabulate details ####################################
	#############################################################################################################	
	'''	

	def download_ALMA_srcMetadata(Src):
		issim=1
		Sno=str(Src[0])
		src=Src[1].replace('\t','')
		print(Sno,' ',src)
		tfNm=Outtabdir+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')
		if Is_SIMBAD_names ==False:
			try:
				simoj=Simbad.query_object(src)
				sid=simoj['MAIN_ID'].data.data[0]
				if 'BYTES' in str(type(sid)).upper():
					sid=sid.decode('UTF-8')
				elif 'STR' in str(type(sid)).upper():
					pass
				else:
					issim=0
					sid=src		
			except:
				issim=0
				sid=src
		'''		
		collect deatls of the star from Good & Bad ALMA obs.
		'''		
		foudflg=0				
		if issim==1:
			dirs=['Good_ones/','Bad_ones_withSIMBADID/']
			allIDs=[]
			locs=[]
			for Dr in dirs:
				print('Analysing from ',Dr,' for star ',src)
				IDar=allID_match[Dr]
				AIDs=[]
				for i in IDar:
					allIDs+=[i.split('@_')[0]]
					AIDs+=[i.split('@_')[0].upper()]
					locs+=[int(i.split('@_')[1])]
				locmatch=dict(zip(allIDs,locs))			
				if src.replace(' ','').upper() in AIDs:
					print('Found ',src,' in ',Dr)
					dets=pickle.load(open(allpf_match[Dr][locmatch[src.replace(' ','')]],'rb'))

					if not os.path.isfile(Outtabdir+Sno+'$_$'+allpf_match[Dr][locmatch[src.replace(' ','')]].split('/')[-1]):
						pickle.dump(dets,open(Outtabdir+Sno+'$_$'+allpf_match[Dr][locmatch[src.replace(' ','')]].split('/')[-1],'wb'))
					else:
						print(Outtabdir+Sno+'$_$'+allpf_match[Dr][locmatch[src.replace(' ','')]].split('/')[-1],' exists!! So not saving metadata')
					if os.path.isdir(tfNm):
						print(tfNm,' exists too!!')
					else:
						os.mkdir(tfNm)
					if not os.path.isdir(tfNm+'/'+tfNm+'_jpg_imgs'):
						os.system('mkdir '+tfNm+'/'+tfNm+'_jpg_imgs')
					if not os.path.isdir(tfNm+'/'+tfNm+'_JVO_imgs'):
						os.system('mkdir '+tfNm+'/'+tfNm+'_JVO_imgs')
					if not os.path.isdir(tfNm+'/'+tfNm+'_JVO_spectra'):
						os.system('mkdir '+tfNm+'/'+tfNm+'_JVO_spectra')

					if os.path.isfile(Outtabdir+Sno+'$_$'+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')+'_allurls.txt'):
						print('All urls are already saved in ',Outtabdir+Sno+'$_$'+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')+'_allurls.txt')
					else:
						frtxtfil='\n'.join(dets[6])
						fil=open(Outtabdir+Sno+'$_$'+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')+'_allurls.txt','a')
						fil.write(frtxtfil)
						fil.close()
					foundflg=1	
					del dets
					break
					
		else:
			if src.replace('/','_slh_').replace(' ','') in Bad_noIDs_nosp:
				tff=allworstpfmatch[src.replace(' ','').replace('/','_slh_')]
				dets=pickle.load(open(tff,'rb'))
				if not os.path.isfile(Outtabdir+Sno+'$_$'+tff.split('/')[-1]):
					pickle.dump(dets,open(Outtabdir+Sno+'$_$'+tff.split('/')[-1],'wb'))
				if not os.path.isfile(Outtabdir+Sno+'$_$'+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')+'_allurls.txt'):
					frtxtfil='\n'.join(dets[-1])
					fil=open(Outtabdir+Sno+'$_$'+src.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_')+'_allurls.txt','a')
					fil.write(frtxtfil)
					fil.close()
				os.system('mkdir '+tfNm)
				os.system('mkdir '+tfNm+'/'+tfNm+'_jpg_imgs')
				foundflg=1								
				del dets				
			if foundflg==0:
				print('No ALMA entry for ',src)
				pickle.dump([src,'--',['--'],['--'],['--'],['--'],['--'],['--'],['--']],open(Outtabdir+Sno+'$_$'+src.replace('/','_slh_').replace(' ','_').replace(';','_smcln_')+'.p','wb'))
		return
		
	'''
	#############################################################################################################	
	#### Download ALMA fits images & generate jpgs of the sources required #####################################
	#############################################################################################################	
	'''	
		
	def Down_ALMAdatafits(tup):
		tffil=tup[0].split('/')[-1]
		tdest=tup[1]
		print('\n__________________________\n')
		print('Downloading .....')
		print('Star: ',tdest)
		print('ID: ',tffil)
		print('\n__________________________\n')
		#jptag=tdest+'_jpg_imgs/'
		jvo_img_loc=tdest+'_JVO_imgs/'
		jvo_spec_loc=tdest+'_JVO_spectra/'
		jptag=''
		tdest+='/'
		tdest=Outtabdir+tdest
		jpgdest=tdest+jptag
		jil=tdest+jvo_img_loc
		jsl=tdest+jvo_spec_loc		
		beg=time.time()
		sts='Fits not requested to download!!'
		jvsts=''
		tg=''
		anyIfile=glob.glob(jpgdest+'*'+tffil+'.jpg')
		anyFfile=glob.glob(tdest+tffil+'.fits')
		cubetag=''
		if Get_only_Quicklookimg==False:
			sts='Fits download failed!!'	
			filo=open(tdest+'Bad_IDs.txt','a')				
			if len(anyIfile)==0:
				if not len(anyFfile)==0:
					os.system('rm -rf '+anyFfile[0])
				try:	
					#Downloading using urllib request
					urllib.request.urlretrieve(tup[0],tdest+tffil+'.fits') 
					Img_data=fits.getdata(tdest+tffil+'.fits')[0][0]					
					thd=fits.getheader(tdest+tffil+'.fits')			
					if thd[NAXIS3]>1:
						cubetag='Cube_'
						os.system('mv '+tdest+tffil+'.fits '+tdest+cubetag+tffil+'.fits')
					obj=thd['OBJECT']
					fq=np.round(thd['RESTFRQ']*10**-9,2)
			
					if np.nanmax(Img_data)<=0:
						filo.write('\n'+tffil+'\t Max flux in image is <=0\t'+str(fq))
						print('Bad data for ',tup[1],' since max value of the image data <0')
						tg='Bad_'
					else:
						sts='Download succeeded!!!'	
					plt.figure(figsize=(1200/my_dpi, 800/my_dpi), dpi=my_dpi)	
					plt.imshow(Img_data,aspect='auto',origin='lower',cmap='jet')
					imMax=np.nanmax(Img_data)
					Isize=Img_data.shape[0]
					imRMS=np.nanstd(Img_data[5:Isize-10,5:int(Isize/4)])
					imDR=np.round(imMax/imRMS,3)
					plt.title(obj+' Max: '+str(np.round(imMax,5))+' DR: '+str(imDR)+' RMS err: '+str(np.round(imRMS,8)),fontsize=16,fontweight='bold')
					plt.figtext(0.6,0.9,'Frq: '+str(fq)+' GHz',color='r',fontsize=16,fontweight='bold')
					plt.colorbar()
					plt.tight_layout()
					plt.savefig(jpgdest+tg+cubetag+tffil+'.jpg')
					plt.close()
					print('Jpg image saved for ',tup[1])
					del Img_data
					del thd
				except:
					print('Bad data for ',tup[1])
					filo.write('\n'+tffil+'\tFailed Download\t--')
				filo.close()
			else:
				print(anyIfile[0],' already exists so moving to next star ID.')
				filo.close()
				sts='Fits Data was already downloaded'
		try:
			urllib.request.urlretrieve(jvo_imgurl+tffil,jil+tffil+'.jpg')
			jvsts='JVO image downloaded!!\t'
		except:
			print('JVO image download failed for star ',tup[1])
		try:
			urllib.request.urlretrieve(jvo_specurl+tffil,jsl+tffil+'.jpg')
			jvsts+='JVO Spectrum obtained!!'
		except:
			print('JVO spectra download failed for star ',tup[1])					
		endt=(time.time()-beg)/60.
		print('\n__________________________\n')
		print('Completed analysis .....')
		print('Star: ',tup[1])
		print('ID: ',tffil.replace('.fits',''))
		print('FITS Status: ', sts)
		print('JVO data status: ',jvsts)
		print('Time taken: ',endt,' min.')
		print('\n__________________________\n')
		return		
	########################################################################################		

	if Search_stars_meta==True and just_get_SearchStarimgs==False:
		#os.system('mkdir Temp_files/')
		search_stars_list=[[i,search_stars[i]] for i in range(len(search_stars))]		
		p=Pool(Nproc)
		os.nice(Nice)
		p.map(download_ALMA_srcMetadata,search_stars_list)
		os.chdir(Outtabdir)
		fil_al=sorted(glob.glob('*.p'))
		print(fil_al)
		fil=open(Out_table_name+'.tsv','a')
		fil.write('SNo\tSIMBAD ID\tSource names\tBands\tComments\tALMA_ID\tObs_datteime\tFits file count\tGood FITS files\tGood Extns')
		tol_exp=len(search_stars)
		snos_pst=[]
		for fl in fil_al:
			dets=pickle.load(open(fl,'rb'))
			tsn=fl.split('$_$')[0]
			snos_pst+=[int(tsn)]
			Gc=np.where(np.array(dets[7])!='t')[0]
			goodfils=np.array(dets[8])[Gc]
			Glocs=[]
			goodextns=[]
			for j in range(len(Gc)):
				Fflg=0
				if '.image.' in goodfils[j] or '.pbcor.' in goodfils[j]:
					Glocs+=[Gc[j]]
					Fflg=1
				elif '.flux' in goodfils[j] or '.fits' not in goodfils[j] or '.psf' in goodfils[j] or '.sumwt' in goodfils[j] or '.model' in goodfils[j] or '.residual' in goodfils[j] or '.mask.' in goodfils[j] or '.pb.' in goodfils[j]:
					pass
				else:
					Glocs+=[Gc[j]]
					Fflg=1
				if Fflg==1:	
					tgf=dets[8][Gc[j]].split('.')
					if len(tgf)>=4:
						goodextns+=[tgf[-3]+'.'+tgf[-2]+'.fits']
					elif len(tgf)==3:
						goodextns+=[tgf[-2]+'.fits']
					elif len(tgf)==2:
						goodextns+=[tgf[-1]]	
			pnm=fl.replace(' ','_').replace('*','_star_').replace(';','_smcln_').replace('/','_slh_').replace('.p','_goodObs.p')		
			Ngs=str(len(Glocs))
			pickle.dump([Glocs,list(np.array(dets[2])[Glocs]),list(np.array(dets[6])[Glocs])],open(pnm,'wb'))# Good data locations, ALMA IDs, URLs
			
			tfwrt=tsn+'\t'+dets[0]+'\t'+dets[1]+'\t'+';'.join(dets[4])+'\t'+';'.join(dets[5])+'\t'+';'.join(dets[2])+'\t'+';'.join(dets[3])+'\t'+str(len(dets[6]))+'\t'+Ngs+'\t'+';'.join(goodextns)
			fil.write('\n'+tfwrt)
		fil.close()
		bad_snos=list(set(np.arange(tol_exp))-set(snos_pst))
		bad_vals=search_stars[bad_snos]
		filb=open('Unfound_stars.txt','a')
		filb.write('\n'.join(bad_vals))
		filb.close()
		alllists=glob.glob('*$_$*goodObs.p')
		all_dests=[]
		all_urls=[]
		n_urls=[]
		for fi in alllists:
			tdata=pickle.load(open(fi,'rb'))
			all_urls+=tdata[-1]
			fi_d=fi.replace('_goodObs.p','').split('$_$')[1]
			all_dests+=[fi_d]*len(tdata[0])
			if Download_all==True:
				filb=open(fi_d+'/Bad_IDs.txt','w')
				filb.write('ID\tIssue\tFreq.')
				filb.close()

		url2fil=[(all_urls[i],all_dests[i]) for i in range(len(all_urls))]
		print('#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n____________________________________\n')
		print('Total no. of fits files: ',len(all_urls))
		pickle.dump(url2fil,open('All_urls_with_folderNames.p','wb'))
		print('\n______________________________________\n#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n')
		if Download_all==True:
			
			p.map(Down_ALMAdatafits,url2fil)
	if just_get_SearchStarimgs==True:
		p=Pool(Nproc)
		os.nice(Nice)
		os.chdir(Outtabdir)
		alllists=glob.glob('*$_$*goodObs.p')
		for fi in alllists:
			tdest=fi.replace('_goodObs.p','').split('$_$')[1]
			jvo_img_loc=tdest+'_JVO_imgs/'
			jvo_spec_loc=tdest+'_JVO_spectra/'
			jptag=tdest+'_jpg_imgs/'
			tdest+='/'
			tdest=Outtabdir+tdest
			jpgdest=tdest+jptag
			jil=tdest+jvo_img_loc
			jsl=tdest+jvo_spec_loc			
			if not os.path.isdir(jpgdest):
				os.mkdir(jpgdest)
			if not os.path.isdir(jil):
				os.mkdir(jil)
			if not os.path.isdir(jsl):
				os.mkdir(jsl)			
			os.system('mv '+tdest+'*.jpg '+jpgdest)

		try:
			url2fil=pickle.load(open('All_urls_with_folderNames.p','rb'))
		except:
			alllists=glob.glob('*$_$*goodObs.p')
			all_dests=[]
			all_urls=[]
			n_urls=[]
			for fi in alllists:
				tdata=pickle.load(open(fi,'rb'))
				all_urls+=tdata[-1]
				fi_d=fi.replace('_goodObs.p','').split('$_$')[1]
				all_dests+=[fi_d]*len(tdata[0])
				filb=open(fi_d+'/Bad_IDs.txt','w')
				filb.write('ID\tIssue\tFreq.')
				filb.close()

			url2fil=[(all_urls[i],all_dests[i]) for i in range(len(all_urls))]

			print('#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n____________________________________\n')
			print('Total no. of fits files: ',len(all_urls))
			pickle.dump(url2fil,open('All_urls_with_folderNames.p','wb'))
			print('\n______________________________________\n#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n')
		p.map(Down_ALMAdatafits,url2fil)	
#################################################################################################################

################################################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
################################################################################################################################	

'''	
#############################################################
#VIZIER & SIMBAD based recording of stellar properties. Also estimates the unknown parameters using other 
#parameters based on various correlation functions.
##################################################################	
'''
		
def info_search(st):
	print('\n------------------------')
	print('\n************************\n')
	print('Collecting data for star: ',st)
	print('\n************************\n')
	print('\n------------------------\n')

	beg=time.time()
	Keys=[]
	vals=[]
	succ_cats=[] # All successful catalogs
	rad_min=np.max([radius_min,0.5])
	rad_max=np.min([radius_max,6.])
	del_r=(rad_max - rad_min)/5.
	bfil=open(basedir+'/Logs/'+MY_catalog+'_No_Vizieroutput.txt','a')
	if del_r>.5:
		rads=np.linspace(rad_min,rad_max,5)
	else:
		rads=np.arange(rad_min,rad_max,0.5)	
	for radius in rads:
		try:	
			obj=Vizier.query_object(st,radius=radius*u.arcsec) #Query Vizier
		except:
			try:
				obj=Vizier.query_object(st,radius=3*u.arcsec)
				print('\nFailed Vizier query with usual route. So used radius=3 arcsec and succeeded.')				
			except:	
				print('Vizier Query error for star,', st ,'..')
				bfil.write(st+'\n')
				bfil.close()
				#return -1
		if len(obj.keys())>0:
			break
	if len(obj.keys())<0:
		obj=Vizier.query_object(st,radius=5*u.arcsec)
	Tcats=obj.keys()	#Loading all Catalogs with the star info
	imp=list(Catalog_pref.values())	#List of user specified important catalogs
	# Checking if user specified catalogs have the star info. If yes placing those catalogs on top of the search list.
	catalogs=list(set(imp).intersection(set(Tcats)))
	if look_across_cats==True:
		catalogs+=list(set(Tcats)-set(imp))
	try:	
		simoj=Simbad.query_object(st) # SIMBAD data collection	
	except:
		print('SIMBAD Query error for star,', st ,'..')
		bfil.write(st+'\n')
		bfil.close()
		return -1
#############################################################################
#Searching VIZIER for the value of properties specified by user.	
##############################################################################
	bfil.close()
	bv_flg=0
	bad_keys=[]
	calcflg=0
	for prk,prv1 in Props.items():
		print('\n*******************************\n','Checking for ',prk,' for star ',st,'\n*******************************\n')
		prv=prv1[0]
		prv_typ=prv1[1]
		prv_ucds=prv1[2]
		stflg=0 # Flag to notify the dtype expected in string
		numflg=0
		if 'str' in key_typ[prk]:
			stflg=1
		if 'int' in key_typ[prk] or 'float' in key_typ[prk]:
			numflg=1	
		#print('Finding ',prk,' with keywords: ',prv)
		#print('Value expected: ',prv_typ,'\nUCD: ',prv_ucds)
		Keys+=[prk]
		Cat_entry='--'
		Cat='--'
		Cat_val='--'
		Cat_unit='--'
		flg=0
		for ct in catalogs:
			ctitems=obj[ct].keys()
			#print('Checking Catalog: ',ct,' for parameter, ',prk,' for star ',st)
			for j in ctitems:
				jt=''.join(filter(str.isalnum,j)).upper()
				tatt=[''.join(filter(str.isalnum,att)).upper() for att in prv ]
				#print('Item: ',j)
				#print('jt: ',jt,'\ntatt: ',tatt)
				if jt in tatt:
					### Standardising the units
					Cat_val=obj[ct][j].data[0]
					### If 'ucd' key isnt there in the metadata then its most likely unbelievable or wrong info.
					if 'ucd' in obj[ct][j].meta.keys():
						cat_ucd=obj[ct][j].meta['ucd'].replace(' ','')
					else:
						continue
					###########################################################################################
										
					### If string or number is the expected keyword value, then check if the expected type match the actual.
					if stflg==1 and 'STR' not in str(type(Cat_val)).upper():
						#print(type(Cat_val),' mismatch str. So trying other catalog items.')
						continue
					if numflg==1:	
						if 'FLOAT' in str(type(Cat_val)).upper() or 'INT' in str(type(Cat_val)).upper():
							pass
						else:							
							#print(type(Cat_val),' mismatch number. So trying other catalog items.')	
							continue
					##############################################################################################		
					
					#print('Catalog UCD: ',cat_ucd,'\n Value of Param: ',Cat_val,'\nType of param: ',type(Cat_val),'\nAllowed UCDs: ',prv_ucds)
					#print(cat_ucd in prv_ucds)
					#print('Matching with expected dtype: ',prv_typ(Cat_val))
					tcat_val=''.join(filter(str.isalnum,str(Cat_val)))
					#print('Found ',prk,' in ',ct,' as ',j)	
					# Checking if the property has a value entered, if the property name is one of the different user specified
					# property identifiers, if the datatype of the value is correct and if property UCD also match with user specified value.
					if str(Cat_val)!='None' and tcat_val not in tatt and prv_typ(Cat_val) and cat_ucd in prv_ucds:
						Cat_entry=j
						Cat=ct
						if prk=='M':
							if "solMass" in str(obj[ct][j].unit) or str(obj[ct][j].unit)==u.Msun or "Msun" in str(obj[ct][j].unit):
								Cat_val=obj[ct][j].data[0]
								Cat_unit='Msun'
							elif not obj[ct][j].unit == u.Msun:
								dtt=str(type(obj[ct][j].unit)).upper()
								if 'NONE' in dtt or 'UNRECOGNIZED' in dtt:
									print('Bad_units for',prk,' so next itreation..!!')
									continue	
								else:
									tvl=obj[ct][j].data[0]*obj[ct][j].unit
									try:
										Cat_val=tvl.to(u.Msun).value
										Cat_unit='Msun'
										print('Unit of mass changed to Msun')				
									except:
										Cat_val=obj[ct][j].data[0]
										Cat_unit=str(obj[ct][j].unit)
										
							elif obj[ct][j].data[0]<0:
								print('Found -ve value for ',prk,' with correct units. So taking 10^x..!!')							
								Cat_val=np.round(10**obj[ct][j].data[0],val_prec)
								Cat_unit='Msun'
						if prk=='R':
							if "solRad" in str(obj[ct][j].unit) or str(obj[ct][j].unit)==u.Msun or "Rsun" in str(obj[ct][j].unit):
								Cat_val=obj[ct][j].data[0]
								Cat_unit='Rsun'
							elif not obj[ct][j].unit == u.Rsun:
								dtt=str(type(obj[ct][j].unit)).upper()
								if 'NONE' in dtt or 'UNRECOGNIZED' in dtt:
									print('Bad_units for',prk,' so next itreation..!!')
									continue	
								else:
									tvl=obj[ct][j].data[0]*obj[ct][j].unit
									try:
										Cat_val=tvl.to(u.Rsun).value
										Cat_unit='Rsun'
										print('Unit of radius changed to Rsun')				
									except:
										Cat_val=obj[ct][j].data[0]
										Cat_unit=str(obj[ct][j].unit)
										
							elif obj[ct][j].data[0]<0:
								print('Found -ve value for ',prk,' with correct units. So taking 10^x..!!')							
								Cat_val=np.round(10**obj[ct][j].data[0],val_prec)
								Cat_unit='Rsun'


						if prk=='Lbol':
							if "solLum" in str(obj[ct][j].unit) or "Lsun" in str(obj[ct][j].unit) or obj[ct][j].unit==u.Lsun:
								Cat_val=obj[ct][j].data[0]
								Cat_unit='Lsun'
							elif not obj[ct][j].unit == u.Lsun:
								dtt=str(type(obj[ct][j].unit)).upper()
								if 'NONE' in dtt or 'UNRECOGNIZED' in dtt:
									print('Bad_units for',prk,' so next itreation..!!')
									continue	
								else:
									tvl=obj[ct][j].data[0]*obj[ct][j].unit
									try:
										Cat_val=tvl.to(u.Lsun).value
										Cat_unit='Lsun'
										print('Changed ',prk,' to Lsun')
									except:
										Cat_val=obj[ct][j].data[0]
										Cat_unit=str(obj[ct][j].unit)
								
							elif obj[ct][j].data[0]<0:
								print('Found -ve value for ',prk,' with correct units. So taking 10^x..!!')
								Cat_val=np.round(10**obj[ct][j].data[0],val_prec)
								Cat_unit='Lsun'
						if prk=='Distance':
							#print('Found distance')
							if obj[ct][j].unit == u.pc or "pc" in str(obj[ct][j].unit) or obj[ct][j].unit == u.kpc or "kpc" in str(obj[ct][j].unit):
								Cat_val=obj[ct][j].data[0]
								if obj[ct][j].unit == u.kpc or "kpc" in str(obj[ct][j].unit):
									Cat_val=obj[ct][j].data[0]*1000
								Cat_unit='pc'								
							elif not obj[ct][j].unit == u.pc:
								dtt=str(type(obj[ct][j].unit)).upper()
								if 'NONE' in dtt or 'UNRECOGNIZED' in dtt:
									print('Bad_units for',prk,' so next itreation..!!')
									continue	
								else:
									tvl=obj[ct][j].data[0]*obj[ct][j].unit
									try:
										Cat_val=tvl.to(u.pc).value
										Cat_unit='pc'
										print('Distance unit changed to pc')
									except:
										Cat_val=obj[ct][j].data[0]
										Cat_unit=str(obj[ct][j].unit)

							#print('unit: ',Cat_unit)	
						if prk=='Age':
							#print('Found age in',ct)
							if "Myr" in str(obj[ct][j].unit) or obj[ct][j].unit == u.Myr or "Gyr" in str(obj[ct][j].unit):
								Cat_val=obj[ct][j].data[0]
								if "Gyr" in str(obj[ct][j].unit):
									Cat_val=obj[ct][j].data[0]*1000
								Cat_unit='Myr'
							elif not obj[ct][j].unit == u.Myr:
								dtt=str(type(obj[ct][j].unit)).upper()
								if 'NONE' in dtt or 'UNRECOGNIZED' in dtt:
									print('Bad_units for',prk,' so next itreation..!!')
									continue	
								else:
									print('Wrong unit, ',str(obj[ct][j].unit))
									tvl=obj[ct][j].data[0]*obj[ct][j].unit
									try:
										Cat_val=tvl.to(u.Myr).value
										Cat_unit='Myr'
									except:
										Cat_val=obj[ct][j].data[0]
										Cat_unit=str(obj[ct][j].unit)
						if prk=='logRx':
							if obj[ct][j].data[0]>0:
								Cat_val=np.round(np.log10(obj[ct][j].data[0]),val_prec)
						if prk=='logLx':
							if jt in ['LX','LXRAY']:
								if obj[ct][j].data[0]>100000:
									Cat_val=np.round(np.log10(obj[ct][j].data[0]),val_prec)			
								elif obj[ct][j].data[0]>0:
									Cat_val=obj[ct][j].data[0]
								else:
									continue	
						print('\n------------------\nPassed the gate..\n')
						flg=1
						if Cat not in succ_cats:
							print('Adding Catalog ',Cat,' to dictionary.')
							succ_cats+=[Cat]						
						#Cat_val=str(Cat_val)
						print('Recorded ',prk,' value for ',st,' : ',Cat_val,' in ',ct)
						#print('Found ',tatt,' for ',st,' in ',j)
						break
					else:
						Cat_val='--'
						Cat_unit='--'
						flg=0
				if flg==1:
					break								 
			if flg==1:
				break
		if flg==1:
			if prk=='B-V':
				Cat_unit='mag'
				if Cat_val<=1.5 and Cat_val>=0.5:
					bv_flg=1	
			if prk=='V':
				Cat_unit='mag'
			if Cat_unit=='--':
				Cat_unit=str(obj[ct][Cat_entry].unit)
				if '[' in Cat_unit:
					Cat_unit='--'							
			vals+=[[Cat_val,Cat_unit,Cat]]
			print('Values for ',st,' : ',[Cat_val,Cat_unit,Cat])
			print('\n------------------\n')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Using SIMBAD to get Otype and unknowns
####################################################################################

		if flg==0 and prk in ['B-V','Distance','RA','DEC','R','SpType']:
			#print('Checking SIMBAD for ',prk)
			if prk=='B-V':
				BVmag=simoj['FLUX_B'].data.data[0]-simoj['FLUX_V'].data.data[0]
				BVstr=str(abs(BVmag)).replace('.','')
				if BVstr.isnumeric():
					flg=1
					vals+=[[BVmag,'mag','SIMBAD']]
					print('Found ',prk,' in SIMBAD for ',st)
					if BVmag<=1.5 and BVmag>=0.5:
						bv_flg=1
			if prk == 'Distance':
				tdis=simoj['Distance_distance'].data.data[0]
				tdisstr=str(tdis).replace('.','')
				if tdisstr.isnumeric():
					flg=1
					vals+=[[tdis,simoj['Distance_unit'].data.data[0],'SIMBAD']]							
			if prk in ['RA','DEC']:
				tdis=simoj[prk].data.data[0]
				if prk == 'RA':
					flg=1
					vals+=[[tdis,'hms','SIMBAD']]							
				else:
					flg=1
					vals+=[[tdis,'dms','SIMBAD']]
			if prk == 'R':
				tdis=simoj['Diameter_diameter'].data.data[0]
				tdisstr=str(tdis).replace('.','')
				if tdisstr.isnumeric():
					flg=1
					vals+=[[tdis/2.,simoj['Diameter_unit'].data.data[0],'SIMBAD']]
			if prk=='SpType':
				tdis=simoj['SP_TYPE'].data[0]
				if 'BYTES' in str(type(tdis)).upper():
					if tdis.isalnum():
						tdis=tdis.decode("utf-8")
						flg=1
						vals+=[[tdis,'--','SIMBAD']] 
				elif 'STR' in str(type(tdis)).upper():
					if tdis.isalnum():
						flg=1
						vals+=[[tdis,'--','SIMBAD']]
		if not flg==1:
			vals+=[[Cat,Cat,Cat]]
		#print('Keys:',Keys,'\nLength: ',len(Keys),'\nVals: ',vals,'\nLength: ',len(vals))	
	Keys+=['Other OTypes']
	tmpt=simoj['OTYPES'].data.data[0]
	if 'BYTES' in str(type(tmpt)).upper():
		tmpt=tmpt.decode("utf-8")
	vals+=[[tmpt,'--','SIMBAD']]
	succ_cats+=['SIMBAD']
	Keys+=['SIMBAD OType']
	tmpt=simoj['OTYPE'].data.data[0]
	if 'BYTES' in str(type(tmpt)).upper():
		tmpt=tmpt.decode("utf-8")
	vals+=[[tmpt,'--','SIMBAD']]
	Keys+=['SIMBAD MainID']
	tmpt=simoj['MAIN_ID'].data.data[0]
	if 'BYTES' in str(type(tmpt)).upper():
		tmpt=tmpt.decode("utf-8")
	vals+=[[tmpt,'--','SIMBAD']]
	
	print('Star ',st,' OType: ',simoj['OTYPE'].data.data[0],'\nMAIN ID: ',simoj['MAIN_ID'].data.data[0])
			#print('Values: ',[Cat,Cat,Cat])
	# Making dictionary of Property:[Value,Unit,Catalog]
	print('Keys:',Keys,'\nLength: ',len(Keys),'\nVals: ',vals,'\nLength: ',len(vals))
	Dct=dict(zip(Keys,vals))

###################################################	
# Calculate Unknown Parameters ####################
###################################################	

	if Calc_params==True:
		print('Calculating params for ',st)
		if 'Mv' in Keys:
			vmagT=simoj['FLUX_V'].data.data[0]
			vstr=str(abs(vmagT)).replace('.','')
			if not Dct['Distance'][0]=='--' and Dct['Distance'][1]=='pc' and vstr.isnumeric():
				dpc=Dct['Distance'][0]				
				Dct['Mv']=[np.round(vmagT+5*(1-np.log10(dpc)),val_prec),'mag','Calculated']
				calcflg=1
				print('Mv found for ',st,' : ',Dct['Mv'][0])
		#######################################################		
		#if M or Teff or Lbol is there then use it to find B-V.
		#######################################################
		BVmagt=[-99,0,0]
		if not Dct['M'][0]=='--':
				if bv_flg==0:
					#print(Dct['M'][0])
					BVmagt=[np.round((0.28-np.log10(Dct['M'][0]))/0.42,val_prec),'mag','Calculated']
					calcflg=1
					bv_flg=2
				if not Dct['Lbol'][0]=='--':
					calcflg=1
					Dct['Lbol']=[np.round(MtoL(Dct['M'][0]),val_prec),'Lsun','Calculated']
		elif not Dct['Lbol'][0]=='--':
			Dct['M']=[np.round(LtoM(Dct['Lbol'][0]),2),'Msun','Calculated']
			calcflg=1
			if bv_flg==0:
				calcflg=1
				BVmagt=[np.round((0.28-np.log10(Dct['M'][0]))/0.42,val_prec),'mag','Calculated']
				bv_flg=2
		elif not Dct['Teff'][0]=='--':
			if bv_flg==0:
				calcflg=1		
				BVmagt=[np.round((3.908-np.log10(Dct['Teff'][0]))/0.234,val_prec),'mag','Calculated']
				bv_flg=2
		if bv_flg==2 and BVmagt[0]<=1.5 and BVmagt[0]>=0.5:
			bv_flg=1
			calcflg=1
			Dct['B-V']=BVmagt	
		del BVmagt		
		#######################################################################	
		if bv_flg==1:
			all_ryt=False
			if not Dct['logRHK'][0]=='--':
				calc_vals=Spec_Prop(Dct['B-V'][0],rhkp=10**Dct['logRHK'][0])
				Dct['S']=[np.round(calc_vals.S,val_prec),'--','Calculated']
				calcflg=1
				all_ryt=True
			elif not Dct['S'][0]=='--':
				calc_vals=Spec_Prop(Dct['B-V'][0],S=Dct['S'][0])
				Dct['logRHK']=[np.round(np.log10(calc_vals.Rhkp),val_prec),'--','Calculated']
				calcflg=1
				all_ryt=True
			if all_ryt==True:
				if Dct['M'][0]=='--':
					Dct['M']=[np.round(calc_vals.M,val_prec),'Msun','Calculated']
					calcflg=1
				if Dct['Teff'][0]=='--':
					calcflg=1
					Dct['Teff']=[np.round(calc_vals.M,val_prec),'K','Calculated']
				if 'Period' in Keys and Dct['Period'][0]=='--':
					Dct['Period']=[np.round(calc_vals.period,val_prec),'d','Calculated']
					calcflg=1
				if 'Ro' in Keys:
					if Dct['Ro'][0]=='--':
						Dct['Ro']=[np.round(calc_vals.R0,val_prec),'--','Calculated']
						calcflg=1
				if 'logRo' in Keys:
					if Dct['logRo'][0]=='--':
						Dct['logRo']=[np.round(np.log10(calc_vals.R0),val_prec),'--','Calculated']
						calcflg=1
				if 'logRx' in Keys:
					if Dct['logRx'][0]=='--':
						Dct['logRx']=[np.round(np.log10(calc_vals.Rx),val_prec),'--','Calculated']
						calcflg=1
				if 'B_l' in Keys:
					if Dct['B_l'][0]=='--':
						Dct['B_l']=[np.round(calc_vals.B_surf,val_prec),'G','Calculated']
						calcflg=1
			else:
		###########################
		#Find M Teff Lbol##########
		###########################
				if Dct['M'][0]=='--':
					Dct['M']=[np.round(10**(0.28-0.42*Dct['B-V'][0]),val_prec),'Msun','Calculated']
					calcflg=1
				if Dct['Teff'][0]=='--':
					Dct['Teff']=[np.round(10**(3.908-0.234*Dct['B-V'][0]),val_prec),'K','Calculated']
					calcflg=1
				if Dct['Lbol'][0]=='--':
					Dct['Lbol']=[np.round(MtoL(Dct['M'][0]),val_prec),'Lsun','Calculated']		
					calcflg=1										
		elif Dct['logRHK'][0] is float:
			y=5+Dct['logRHK'][0]
			if 'logRo' in Keys:
				Dct['logRo']=[np.round(0.324-0.4*y-0.283*y**2-1.325*y**3,val_prec),'dex','Calculated']
				calcflg=1
			if 'Ro' in Keys:
				Dct['Ro']=[np.round(10**(0.324-0.4*y-0.283*y**2-1.325*y**3),val_prec),'--','Calculated']
				calcflg=1
		if Dct['Lbol'][0] is float and Dct['Teff'][0] is float and Dct['R'][0]=='--':
			Ls= 3.828*10**26
			Sbc=5.67*10**-8
			Rsn=695500000

			tLb=Dct['Lbol'][0]*Ls
			tTf=Dct['Teff'][0]
			doit=False
			if 'K' not in Dct['Teff'][1]:
				if 'Tsun' in Dct['Teff'][1]:
					Dct['Teff'][0]=Dct['Teff'][0]*5778
					Dct['Teff'][1]='K'
					doit=True
			else:
				doit=True
			if doit==True:
				tRst=np.sqrt(tLb/(Sbc*tTf**4*4*np.pi))
				Dct['R']=[np.round(tRst/Rsn,val_prec),'Rsun','Calculated']
				calcflg=1	
	if calcflg==1:
		succ_cats+=['Calculated']		
				
#########################################################################################					

	fname=st.replace(' ','_')+'.p'
	#Record the name of star for which data is stored as .p file
	fil=open(basedir+'/Logs/'+MY_catalog+'_Stars_done.txt','a')
	if fil.tell()>0:
		fil.write('\n'+st)
	else:
		fil.write(st)	
	fil.close()
	pickle.dump(succ_cats,open(pfildir_cats+st.replace(' ','_')+'_Catalogs_used.p','wb'))
	pickle.dump(Dct,open(pfildir+fname,'wb'))	
	endt=(time.time()-beg)/60.
	print('Result for star, ',st,'\n',Dct)
	print('Done in',endt,' min.')
	print('Catalogs: ',succ_cats)
	del Keys
	del vals
	del Dct
	return flg

'''
##########################################################################################################
### Tabulating the stellar details collected/determined from/using Vizier and SIMBAD #####################
##########################################################################################################				
'''	
def make_table():
	#pfildir=Table_dir+MY_catalog+'_Stardata/'
	os.chdir(pfildir_cats)
	str_cats=glob.glob('*.p')
	#strs=glob.glob('*.p')
	allstars=[i.replace('_Catalogs_used.p','') for i in str_cats]
	done_with=np.array([])
	calflag=0
	if update_bib==True:
		Super_num,Super_cat=np.loadtxt(Table_dir+bib_txtfile,dtype=str,unpack=1,usecols=(0,1),skiprows=1,delimiter='\t\t')
		Super_num=list(Super_num)
		Super_cat=list(Super_cat)
		bib_file=Table_dir+bib_txtfile
		print('Updating ....',bib_file)
		done_with=np.loadtxt(Table_dir+MY_catalog+'.tsv',dtype=str,unpack=1,usecols=(0),skiprows=1,delimiter='\t')
		print('Stars done with: ',done_with,'\nNumber of stars done recording: ',len(done_with))
		if '*' not in Super_num:
			sno=len(Super_num)+1
		else:
			print('Removing Calculated..')
			sno=len(Super_num)
			Super_num.remove('*')
			Super_cat.remove('Calculated')
			calflag=1	
		print('Catalog count : ',sno)			            
	else:
		Super_cat=[]
		Super_num=[]
		sno=1
		bib_file=Table_dir+MY_catalog+'_References.txt'
	print('Super_cat: ',Super_cat,'Super_num: ',Super_num)
	print('\n****************\nStar count: ',len(allstars),'\nTotal stars Recorded: ',len(done_with))
	rem_strs=list(set(allstars)-set(done_with))
	print('Stars to be recorded: ',rem_strs,'\n Total count: ',len(rem_strs),'\n*********************\n')

	for sc in rem_strs:
		print('Recording catalogs with Star: ',sc)
		tsc=pickle.load(open(sc+'_Catalogs_used.p','rb'))
		for nm in tsc:
			if nm not in Super_cat:
				print('Adding catalog ..',nm)
				Super_cat+=[nm]
				Super_num+=[sno]
				sno+=1
	if 'Calculated' in Super_cat or calflag==1:
		Super_cat.remove('Calculated')
		Super_cat+=['Calculated']
		Super_num=Super_num[:-1]
		Super_num+=['*']			
	Cat_dict=dict(zip(Super_cat,Super_num))
	print('Catalog dictionary: ',Cat_dict)
	os.chdir(Table_dir)
	
	if delete_old == True:
		os.system('rm -rf '+MY_catalog+'.tsv')
	fil=open(MY_catalog+'.tsv','a')
	############################################################
	##### Entering the column titles ###########################
	Props.pop('SpType')
	if is_new_catalog==True:
		fil.write('Name')
		fil.write('\tSIMBAD MAIN ID\tSpType\tSIMBAD OTYPE\tOther Otypes')			
		for k in list(Props.keys()):
			fil.write('\t'+k+'\tUnit')
	fil.write('\n')
	#############################################################################
	## Updating/Making the catalog using stellar data ###########################
	#############################################################################
	for st in rem_strs:
		print('Cataloging Star: ',st)
		Data=pickle.load(open(pfildir+st+'.p','rb'))
		fil.write(st.replace('.p',''))

		dty=Data['SIMBAD MainID'][0]
		if 'BYTES' in str(type(dty)).upper():
			dty=dty.decode('utf-8')		
		fil.write('\t'+dty+' ['+str(Cat_dict[Data['SIMBAD MainID'][-1]])+']')

		if not Data['SpType'][-1]=='--':
			fil.write('\t'+str(Data['SpType'][0])+' ['+str(Cat_dict[Data['SpType'][-1]])+']')		
		else:
			fil.write('\t'+str(Data['SpType'][0]))
			
		dty=Data['SIMBAD OType'][0]
		if 'BYTES' in str(type(dty)).upper():
			dty=dty.decode('utf-8')		
		fil.write('\t'+dty+' ['+str(Cat_dict[Data['SIMBAD OType'][-1]])+']')

		dty=Data['Other OTypes'][0]
		if 'BYTES' in str(type(dty)).upper():
			dty=dty.decode('utf-8')		
		fil.write('\t'+dty+' ['+str(Cat_dict[Data['Other OTypes'][-1]])+']')

		for k in list(Props.keys()):
			#print(Data[k][-1])
			if not Data[k][-1]=='--':
				fil.write('\t'+str(Data[k][0])+' ['+str(Cat_dict[Data[k][-1]])+']\t'+str(Data[k][1]))
			else:
				fil.write('\t'+str(Data[k][0])+'\t'+str(Data[k][1]))					
		fil.write('\n')
		del Data
		###############################################################################################
	fil.close()
	if update_bib==True:
		os.system('rm -rf '+bib_file)
	fil=open(bib_file,'a')
	fil.write('#S.No\t\tReference Table')
	for i,j in Cat_dict.items():
		fil.write('\n'+str(j)+'\t\t'+str(i))
	fil.close()	
	os.chdir(basedir)
	return
'''
##########################################################################################################
### Dowloading useful Stellar catalogs from Vizier using publication/catalog table name ##################
##########################################################################################################				
'''	

def Catalog_downloader():
	catnames,viz=np.loadtxt(INPUT_files_folder+'/'+Catalog_names_file, dtype=str,unpack=1, usecols=(0,1), delimiter='\t')
	if 'array' not in str(type(catnames)):
		catnames=np.array([catnames])
		viz=np.array([viz])

	for ctn,ctl in zip(catnames,viz):
		tcl=Vizier.find_catalogs(ctl)
		if len(tcl)>1:
			ctl='NONE'
		elif len(tcl)==1:
			print('Tables in the chosen catalog: ')
			cid=list(tcl.keys())[0]
			for i,j in enumerate(tcl[cid].tables):
				print(i,' ',j.name,'\t',j.description)
			chc=int(input('Choose a table serial number: '))
			ctl=tcl[cid].tables[chc].name	

		if 'NONE' in ctl.upper() and '/' not in ctl:
			opts=Vizier.find_catalogs(ctn)
			print('Catalog: ',ctn,'\nChoose one of your choice to download.')
			for i,k in enumerate(opts.keys()):
				print(i,' : ',opts[k])
			chc=int(input('Choose a catalog Serial number: '))
			ctl=list(opts.keys())[chc]
			for i,j in enumerate(opts[ctl].tables):
				print(i,' ',j.name,'\t',j.description,'\n')
			chc=int(input('Choose a table serial number to download: '))
			ctl=opts[ctl].tables[chc].name
			print('Downloading the table.......')
		# Download the identified table
		catalog=Vizier.get_catalogs(ctl)
		catalog[0].write(Catalog_dir+ctn+'.tsv',format='ascii.tab',overwrite=True)
	return

'''
#############################################################################################################	
#### Cataloging ALMA science Observations of user's requirements. It can also filter out redundent sources
#### and accumulate all observations of same source intelligently ###########################################
#############################################################################################################	
'''
	
def ALMA_Obs_catalog(nmlc):
	catalog_name=catalog_names[nmlc]
	Keyword_IDs=Keyword_IDs_list[nmlc]
	Keyword_values=Keyword_values_list[nmlc]

	#############################################################
	############# Searching for ALMA observations ###############
	#############################################################
	
	keydict=dict(zip(Keyword_IDs,Keyword_values))
	res=Alma.query(payload=keydict,public=Public,science=Science_only)
	if len(columns_wanted)>1:
		res=res[columns_wanted]
	if len(res)<1:
		print('No data found for keyword ',Keyword_values)
		return 0
	###################################################################################	
	# Convert every table data type to str ############################################
	###################################################################################
	
	titles=res.keys()
	cloc=0
	for tl in titles:
		nxt_flg=0
		dt=str(type(res[tl][0])).replace('_','').upper()
		print(dt,' found: ',tl)
		if any([True if i in dt else False for i in ['STR','BYTES','FLOAT','INT']]):
			print(tl)
			tlen=np.mean([len(str(res[tl][i])) for i in range(np.min([len(res),5]))])
			if tlen>8:
				res[tl]=res[tl].astype('S8192')
			else:
				res[tl]=res[tl].astype('S4096')
											
		elif 'MASKEDARRAY' in dt:
			print(dt,' found: ',tl,'. Entered MASKEDARRAY section.')
			vals=res[tl][:]
			res[tl].name='badcol'
			cdata=[]
			for i in range(len(vals)):
				try:
					tdt=vals[i].data
					tdt=tdt.astype('<U256')
					#print('tdt:',tdt)
					nxt_flg=2
				except:
					print('Removing ',tl,' from saved table. Unknown data type in masked column.')
					nxt_flg=1
					res.remove_column(tl)
					break
				#print('cdata:',cdata)	
				cdata+=['\t'.join(tdt)]
					
		elif 'ARRAY' in dt:
			print(dt,' found: ',tl,'. Entered Array section.')
			vals=res[tl][:]
			res[tl].name='badcol'
			cdata=[]
			for i in range(len(vals)):
				try:
					tdt=vals[i].astype(str)
					nxt_flg=2
				except:
					print('Removing ',tl,' from saved table. Unknown data type in masked column.')
					res.remove_column(tl)
					nxt_flg=1
					break
				cdata+=['\t'.join(tdt)]
		else:
			print(dt,' found: ',tl,'. Entered bad side.')
			print('Removing ',tl,' from saved table. Unknown data type in masked column.')
			res.remove_column(tl)
			nxt_flag=1
		if nxt_flg==2:			
				res.remove_column('badcol')
				coldat=Column(data=cdata,name=tl)
				res.add_column(coldat,index=cloc)
				res[tl]=res[tl].astype('S2048')						
		cloc+=1
	sources=res['Source name'].data.data
	###########################################################################################
	### Making new column of Frequency ranges and resolution using Frequency Support column####
	###########################################################################################
	
	fqs_loc=-100 # Location of 'Frequency support' keyword in ALMA catalog  
	fqs_loc=np.where(np.array(titles)=='Frequency support')[0][0]
	fqres_loc=np.where(np.array(titles)=='Frequency resolution')[0][0]
	all_fqrange=[]
	all_fqres=[]	
	fs_vals=res['Frequency support'][:].data
	for vl in fs_vals:
		fq_range=''
		fq_res=''	
		tvl=vl.split(' U ')
		for ttl in tvl:
			fq_range+=','+ttl[1:-1].split(',')[0]
			fq_res+=','+ttl[1:-1].split(',')[1]
		fq_range='['+fq_range[1:]+']'
		fq_res='['+fq_res[1:]+']'
		all_fqrange+=[fq_range]
		all_fqres+=[fq_res]

	res.remove_column('Frequency resolution')
	coldat=Column(data=all_fqres,name='Frequency resolution')
	res.add_column(coldat,index=fqres_loc)
	res['Frequency resolution']=res['Frequency resolution'].astype('S8192')	

	coldat=Column(data=all_fqrange,name='Frequency range')
	res.add_column(coldat,index=fqs_loc)
	res['Frequency range']=res['Frequency range'].astype('S8192')	

	##########################################################################################	
	if filter_repeaters==True:

		#############################################################################
		### Removing Alias names of the star#########################################
		#############################################################################

		sources=res['Source name'].data.data
		orig_srcs=len(sources)
		tsr=len(sources)
		Sources=np.array([i.replace('_','').replace(' ','').upper() for i in sources])
		rws=np.array([],dtype='int64')
		for nmlist in known_alias_names:
			locs=[]
			for nmc in nmlist:
				locs+=list(np.where(Sources==nmc)[0])	
			if len(locs)>1:
				rws=np.append(rws,locs[1:])
				for i in locs[1:]:
					for j in range(len(imp_props)):
						res[imp_props[j]][locs[0]]+=';'+res[imp_props[j]][i]
		res.remove_rows(rws)
		sources=res['Source name'].data.data
		print('Removed ',tsr - len(sources),' sources after alias name filtering')

		#############################################################################
		#### Removing redundant names ###############################################		
		#############################################################################
		tsr=len(sources)
		
		tst=[i.replace('_','').replace(' ','').upper() for i in sources]
		tst=np.unique(np.array(tst))
		Sources=np.array([i.replace('_','').replace(' ','').upper() for i in sources])
		rws=np.array([],dtype='int64')

		for Id in tst:
			locs=np.where(Sources==Id)[0]
			if len(locs)>1:
				rws=np.append(rws,locs[1:])
				for i in locs[1:]:
					for j in range(len(imp_props)):
						res[imp_props[j]][locs[0]]+=';'+res[imp_props[j]][i]
		# Remove redundant rows									
		res.remove_rows(rws)
		
		sources=res['Source name'].data.data					
		print('Removed ',tsr - len(sources),' sources after filtering names with different cases and spl characters.')
		tsr=len(sources)
		
	##############################################################
	# Adding SIMBAD IDs to sources as new column #################
	##############################################################
	SIMBAD_srcs=[] # SIMBAD MAIN IDs for each source
	Good_srcs=[]	# Sources with SIMBAD query result and a valid MAIN Id
	Good_locs=[] # Locs/row numbers of good sources in the table
	Bad_srcs=[]	# Sources with no query result.
	No_MainID=[]	# Sources with no valid MAIN ID
	simids=[]

	if SIMBAD_ID_needed==True:
		tlc=-1
		for srt in sources:
			tlc+=1
			try:
				simoj=Simbad.query_object(srt)
				simidv=simoj['MAIN_ID'].data.data[0]
				if 'BYTES' in str(type(simidv)).upper() or 'STR' in str(type(simidv)).upper():
					pass
				else:
					No_MainID+=[srt]
					simids+=['--']
					print(srt,' has no SIMBAD MAIN ID!!! Ignoring the source...Saving into XXX_No_mainid.p log file')
					continue							
				Good_srcs+=[srt]
				Good_locs+=[tlc]
				SIMBAD_srcs+=[simidv]
				simids+=[simidv]
			except:
				srt=srt.upper()
				if '_STAR_' in srt:
					temsrt=srt.replace('_STAR_','* ')
					try:
						simoj=Simbad.query_object(temsrt)
						simidv=simoj['MAIN_ID'].data.data[0]
						if 'BYTES' in str(type(simidv)).upper() or 'STR' in str(type(simidv)).upper():
							pass
						else:
							No_MainID+=[srt]
							simids+=['--']
							print(srt,' has no SIMBAD MAIN ID!!! Ignoring the source...Saving into XXX_No_mainid.p log file')
							continue							
						Good_srcs+=[srt]
						Good_locs+=[tlc]
						SIMBAD_srcs+=[simidv]
						simids+=[simidv]
					except:						
						print('No response from SIMBAD!! Source will be included in catalog. Saved the name in XXX_SIMBAD_noresponse.p file')
						Bad_srcs+=[srt]
						simids+=['?']
				else:
					print('No response from SIMBAD!! Source will be included in catalog. Saved the name in XXX_SIMBAD_noresponse.p file')
					Bad_srcs+=[srt]
					simids+=['?']
							
		coldat=Column(data=simids,name='SIMBAD Main ID')
		res.add_column(coldat,index=3)

	##########################################################################		
	#### Using SIMBAD MAIN IDs to accumulate same source info into 1 Entry####
	##########################################################################

	if filter_repeats_with_SIMBAD_info==True:
		if SIMBAD_ID_needed==False:
			tlc=-1
			for srt in sources:
				tlc+=1
				try:
					simoj=Simbad.query_object(srt)
					simidv=simoj['MAIN_ID'].data.data[0]
					if 'BYTES' in str(type(simidv)).upper() or 'STR' in str(type(simidv)).upper():
						pass
					else:
						No_MainID+=[srt]
						print(srt,' has no SIMBAD MAIN ID!!! Ignoring the source...Saving into XXX_No_mainid.p log file')
						continue							
					Good_srcs+=[srt]
					Good_locs+=[tlc]
					SIMBAD_srcs+=[simidv]
				except:
					print('No response from SIMBAD!! Source will be included in catalog. Saved the name in XXX_SIMBAD_noresponse.p file')
					Bad_srcs+=[srt]
					
		pickle.dump(Bad_srcs,open(basedir+'/Logs/'+catalog_name+'_SIMBAD_noresponse.p','wb'))
		pickle.dump(No_MainID,open(basedir+'/Logs/'+catalog_name+'_No_mainid.p','wb'))
		SIMBAD_srcs=np.array(SIMBAD_srcs).astype(str)
		Good_srcs=np.array(Good_srcs)
		uniqSIM=np.unique(SIMBAD_srcs)
		Good_locs=np.array(Good_locs)
				
		rws=np.array([],dtype='int64')
		locsrep=[Good_locs[np.where(SIMBAD_srcs==i)[0]] for i in uniqSIM]
		for count in range(len(locsrep)):
			if len(locsrep[count])==1:
				continue
			else:
				locs=locsrep[count]
				rws=np.append(rws,locs[1:])
				for i in locs[1:]:
					res['Source name'][locs[0]]+=';'+sources[i]
					for j in range(len(imp_props)):
						res[imp_props[j]][locs[0]]+=';'+res[imp_props[j]][i]
		# Remove redundant rows									
		res.remove_rows(rws)
		sources=res['Source name'].data.data
		print('Removed ',tsr - len(sources),' sources after SIMBAD ID based filtering..')
	print('\n*****************************\nSource count initially: ',orig_srcs,'\nSource count after filtering repeaters : ',len(sources),'\n****************************\n')	
		############# Saving Catalog ##################################################	

	pickle.dump(res,open(ALMA_catalog_dir+'/'+catalog_name+'.p','wb'))			
	res.write(ALMA_catalog_dir+'/'+catalog_name+'.tsv', format='ascii.tab',overwrite=True)
	return 1
	
#################################################################################################################	
'''	
##################################################################################################################
#################### MAIN Function ##############################################################################
#################################################################################################################	
'''
if __name__=='__main__':
	if Option==1:
		if append_existing_logs==False:# Deleting old log files with same name.
			if os.path.isfile(basedir+'/Logs/'+MY_catalog+'_Stars_done.txt'):
				if delete_old_logs==True:		
					os.system('rm -rf '+basedir+'/Logs/'+MY_catalog+'_Stars_done*.txt')
				else:
					lens=len(glob.glob(basedir+'/Logs/'+MY_catalog+'_Stars_done*.txt'))-1
					os.system('mv '+basedir+'/Logs/'+MY_catalog+'_Stars_done.txt '+basedir+'/Logs/'+MY_catalog+'_Stars_done_v'+str(lens)+'.txt')	
		if Only_make_catalog==False:
			st_mainIDs=[]
			Bads=[]
			Goods=[]
			Reps=[]
			rep_mainIDs=[]
			for st in stars:
				try:
					simd=Simbad.query_object(st)
					tmid=simd['MAIN_ID'].data.data[0]
				except:
					Bads+=[st]
					print('Bad star ',st,' removed from stars list.')
					continue
				if tmid not in st_mainIDs:
					Goods+=[st]
					st_mainIDs+=[tmid]
				else:
					rep_mainIDs+=[tmid]
					Reps+=[st]
					print('Repeat entry ',st,' removed.')	
			matchG=dict(zip(st_mainIDs,Goods))
			matchR=dict(zip(Reps,rep_mainIDs))
			print('Repeater stars and Main IDs: ',matchR)
			good_Reps=[matchG[i] for i in rep_mainIDs]
			matchGR=dict(zip(Reps,good_Reps))
			print('Repeater stars and their good stars in list: ',matchGR)
			pickle.dump(matchGR,open(basedir+'/Logs/'+MY_catalog+'_Repeated_stars_with_Good_stars.p','wb'))			
			pickle.dump(matchR,open(basedir+'/Logs/'+MY_catalog+'_Repeated_stars_with_MAINID.p','wb'))
			pickle.dump(matchG,open(basedir+'/Logs/'+MY_catalog+'_Good_Stars_with_MAINID.p','wb'))		
			p=Pool(Nproc)
			os.nice(Nice)
			succ=p.map(info_search,Goods)
		#for st in stars:
		#	succ=info_search(st)
		make_table()
		os.system("echo 'Table is ready!! says ' "+compd+" | mail -s 'Job done' "+email)
		
	if Option==2:
		Catalog_downloader()
		os.system("echo 'Downloaded Catalogs!! says ' "+compd+" | mail -s 'Job done' "+email)
	if Option==3:
		Nproc=np.min([Nproc,5])
		p=Pool(Nproc)
		os.nice(Nice)
		namelocs=np.arange(len(catalog_names))
		succ=p.map(ALMA_Obs_catalog,namelocs)
		os.system("echo 'Downloaded ALMA Observation metadata for the keywords given and made catalogs!! says ' "+compd+" | mail -s 'Job done' "+email)
	if Option==4:
		os.system("echo 'ALMA observations found for the given stars!! says ' "+compd+" | mail -s 'Job done' "+email)	
