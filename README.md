# ChemCheck

Generate exact 2D chemical structures from names or SMILES.

## Run locally
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn main:app --reload
# http://127.0.0.1:8000/ui
```

API
POST /draw (form-data):
name_or_smiles (str)
fmt (png|svg, default png)
show_indices (bool), show_lone_pairs (bool), show_formal_charges (bool)
w (int), h (int)

Deploy (Render)
Push to GitHub
Create Web Service on Render → select repo
Start command from Procfile: uvicorn main:app --host 0.0.0.0 --port $PORT
Visit /ui
