import pandas as pd
import os
import shutil
df=pd.read_csv("/work/lpdi/users/eline/CID/FUN_1Z97/output/af3ternary/merged_all_nbis.csv")
for id in df.ID:
    path="/work/lpdi/users/eline/CID/FUN_1Z97/output/af3ternary/accepted2/"+str(id)
    if not os.path.isdir(path):
        path="/work/lpdi/users/eline/CID/FUN_1Z97/output/af3ternary/accepted1bis/"+str(id)
        if not os.path.isdir(path):
            print("file not found:"+ id)
            continue
    dest_folder = f"/work/lpdi/users/eline/CID/FUN_1Z97/output/af3ternary/accepted3/{id}"
    os.makedirs(dest_folder, exist_ok=True)
    shutil.copytree(path, dest_folder, dirs_exist_ok=True)
    print(f"copied {path} to {dest_folder}")