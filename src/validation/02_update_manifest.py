"""Update public file inventory and transfer-integrity checksums after an edit."""
from pathlib import Path
import subprocess,csv,hashlib
R=Path(__file__).resolve().parents[2]
files=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=R).decode().split('\0')[:-1]
manifest=R/'reports/reproducibility/FILE_MANIFEST.csv';check=R/'reports/reproducibility/CHECKSUMS.sha256'
with manifest.open('w',newline='')as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['path','bytes']);w.writerows((n,(R/n).stat().st_size)for n in files if R/n not in {manifest,check})
check.write_text(''.join(hashlib.sha256((R/n).read_bytes()).hexdigest()+'  '+n+'\n'for n in files if R/n!=check))
print(len(files),'public files; current manifests rebuilt')
