# parse the xml outputs from plip
import pandas as pd
import xml.etree.ElementTree as ET
import glob
import os

#helper function:
def get_n_hbonds(xml_file, ternary=False): # between binder and target if no ligand and between binder and ligand if ligand is in complex file
  tree = ET.parse(xml_file)
  root = tree.getroot()
  countB=0
  for hb in root.iter('hydrogen_bond'):
    #print(hb.attrib, hb.tag)
    if ternary==False: # scoring a binary complex
      ch1=hb.find('reschain').text
      ch2=hb.find('reschain_lig').text
      #print(ch)
      if ch1=="A" and ch2=="B":
        countB+=1
    else: # if ternary
      ch1=hb.find('reschain').text
      if ch1=="B":
        countB+=1
  return(countB)
    
def get_n_salt_bridges(xml_file, ternary=False):
  tree = ET.parse(xml_file)
  root = tree.getroot()
  countB=0
  for hb in root.iter('salt_bridge'):
    #print(hb.attrib, hb.tag)
    if ternary==False: # scoring a binary complex
      ch1=hb.find('reschain').text
      ch2=hb.find('reschain_lig').text
      #print(ch)
      if ch1=="A" and ch2=="B":
        countB+=1
    else: # if ternary
      """ch1=hb.find('reschain').text
      if ch1=="B":
        countB+=1"""
  return(countB) 

def get_n_pi_stacks(xml_file, ternary=False):
  tree = ET.parse(xml_file)
  root = tree.getroot()
  countB=0
  for hb in root.iter('pi_stack'):
    #print(hb.attrib, hb.tag)
    if ternary==False: # scoring a binary complex
      ch1=hb.find('reschain').text
      ch2=hb.find('reschain_lig').text
      #print(ch)
      if ch1=="A" and ch2=="B":
        countB+=1
    else: # if ternary
      """ch1=hb.find('reschain').text
      if ch1=="B":
        countB+=1"""
  return(countB) 

def get_n_hydrophobic_interactions(xml_file, ternary=False):
  tree = ET.parse(xml_file)
  root = tree.getroot()
  countB=0
  for hb in root.iter('hydrophobic_interaction'):
    #print(hb.attrib, hb.tag)
    if ternary==False: # scoring a binary complex
      ch1=hb.find('reschain').text
      ch2=hb.find('reschain_lig').text
      #print(ch)
      if ch1=="A" and ch2=="B":
        countB+=1
    else: # if ternary
      """ch1=hb.find('reschain').text
      if ch1=="B":
        countB+=1"""
  return(countB)

script_start_time = time.time()
###----
# parser
parser = argparse.ArgumentParser()
parser.add_argument("--xml_files", nargs="+", type=str, help="List of plip xml output files ") # list of pdb files 
parser.add_argument("--ternary", action="store_true", help="set to true if plip was exectuted on ternary complexes")
parser.add_argument("--outdir", required=False, type=str, help="output folder for csv file")
parser.add_argument("--prefix", required=False, type=str, help="prefix for csv file")
args = parser.parse_args()
outdir=args.outdir
prefix=args.prefix
ternary=args.ternary
# find interactions
for xml_file in args.xml_files:
  xml_name=os.path.basename(xml_file)
  id=xml_name.replace("_model_report.xml", "")
  
  data={}
  if ternary:
    prefix="BL_"
  else:
    prefix="AB_"
  data[f"{prefix}hbonds"]=get_n_hbonds(xml_file, ternary)
  data[f"{prefix}salt_bridges"]=get_n_salt_bridges(xml_file, ternary)
  data[f"{prefix}pi_stacks"]=get_n_pi_stacks(xml_file, ternary)
  data[f"{prefix}hydrophobic_interactions"]=get_n_hydrophobic_interactions(xml_file, ternary)
  data["id"]=id
  df=pd.DataFrame(data=data, index=[data["id"]])
  csv_path=f"{outdir}/{prefix}_plip_interactions.csv"
  df.to_csv(csv_path, mode="a", index=False, header=not pd.io.common.file_exists(csv_path))


elapsed_time = time.time() - script_start_time
elapsed_text = f"{'%d hours, %d minutes, %d seconds' % (int(elapsed_time // 3600), int((elapsed_time % 3600) // 60), int(elapsed_time % 60))}"
n_binder=len(args.xml_files)
print(f"Finished Plip profiling scoring for {n_binder} complexes. Script execution took: "+elapsed_text)
