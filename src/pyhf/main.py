import pyhf
import pyhf.readxml  # not imported by default!
import json

# read xml file
spec = pyhf.readxml.parse("PDFandDATA_MXs1/my_measurement.xml", "")

# get workspace
ws = pyhf.Workspace(spec)

with open('MXs1.json', 'w') as f:
    json.dumps(ws, f, indent=2)
