# packs driver/questlhsync into a zip the installer can extract (questlhsync/... , deflate)
import os
import sys
import zipfile

src = os.path.join("driver", "questlhsync")
out = os.path.join("build", "installer", "questlhsync-driver.zip")
if not os.path.isdir(src):
    sys.exit("driver folder is missing: " + src)
os.makedirs(os.path.dirname(out), exist_ok=True)
n = 0
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for root, dirs, files in os.walk(src):
        for name in files:
            full = os.path.join(root, name)
            rel = os.path.relpath(full, "driver").replace("\\", "/")
            z.write(full, rel)
            n += 1
if n == 0:
    sys.exit("driver folder is empty")
print("packed", n, "files into", out)
