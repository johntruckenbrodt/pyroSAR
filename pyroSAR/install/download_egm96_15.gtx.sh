#!/usr/bin/env bash
# download EGM96 geoid model to convert heights with GDAL
cd /usr/share/proj
sudo wget https://download.osgeo.org/proj/vdatum/egm96_15/egm96_15.gtx
sudo chmod 644 egm96_15.gtx
