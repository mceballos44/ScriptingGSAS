from pathlib import Path
import sys

BASE_DIR = Path(__file__).resolve().parent
GSAS_PATH = BASE_DIR / "G2"

sys.path.insert(0, str(GSAS_PATH))

from GSASII import GSASIIpath

# Find the correct binary package for this Python + NumPy version
binary_url = GSASIIpath.getGitBinaryLoc()

print("Binary package:")
print(binary_url)

if binary_url is None:
    raise RuntimeError("No compatible GSAS-II binary package was found.")

# Install binaries into the GSAS-II clone
install_location = GSAS_PATH / "GSASII-bin"

GSASIIpath.InstallGitBinary(
    binary_url,
    str(install_location),
    nameByVersion=True
)

print("Finished installing GSAS-II binaries.")