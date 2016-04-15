
from S1_auxil import OSV

x = OSV("/geonfs02_vol1/S1_OSV/AUX_POEORB", "/geonfs02_vol1/S1_OSV/AUX_RESORB")
# print x.maxdate("POE")
x.update()
