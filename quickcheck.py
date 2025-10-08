import sys
import requests

base = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8000"
ok = 0
total = 0

for line in open("tests_smiles.txt"):
    q = line.strip()
    if not q:
        continue
    total += 1
    resp = requests.post(f"{base}/draw", files={"name_or_smiles": (None, q)})
    if resp.ok:
        ok += 1
        print("OK  ", q)
    else:
        print("FAIL", q, resp.text)

print(f"\nSummary: {ok}/{total} passed")
