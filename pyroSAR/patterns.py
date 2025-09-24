###############################################################################
# Reading and Organizing system for SAR images
# Copyright (c) 2016-2023, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
"""
This file contains regular expressions to identify SAR products.
The pattern 'pyrosar' identifies products in pyroSAR's unified naming scheme.
The names of all other expressions correspond to the classes found in pyroSAR.drivers.
"""
pyrosar = r'(?:.*[/\\]|)' \
          r'(?P<outname_base>' \
          r'(?P<sensor>[A-Z0-9]{1,4})_+' \
          r'(?P<acquisition_mode>[A-Z0-9]{1,4})_+' \
          r'(?P<orbit>[AD])_' \
          r'(?P<start>[0-9T]{15})' \
          r'(?:_(?P<extensions>\w*?)|)' \
          r')_*' \
          r'(?:(?P<polarization>[HV]{2})_' \
          r'(?P<proc_steps>[\w-]*)|)' \
          r'(?P<filetype>(?:.tif|.nc|))$'

ceos_ers = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
           r'(?P<processing_stage_flag>[A-Z])' \
           r'(?P<originator_ID>[A-Z\-]{3})' \
           r'(?P<start_day>[0-9]{8})_' \
           r'(?P<start_time>[0-9]{6})_' \
           r'(?P<duration>[0-9]{8})' \
           r'(?P<phase>[0-9A-Z]{1})' \
           r'(?P<cycle>[0-9]{3})_' \
           r'(?P<relative_orbit>[0-9]{5})_' \
           r'(?P<absolute_orbit>[0-9]{5})_' \
           r'(?P<counter>[0-9]{4,})\.' \
           r'(?P<satellite_ID>[EN][12])' \
           r'(?P<extension>(?:\.zip|\.tar\.gz|\.PS|))$'

ceos_psr1 = r'^LED-ALPSR' \
            r'(?P<sub>P|S)' \
            r'(?P<orbit>[0-9]{5})' \
            r'(?P<frame>[0-9]{4})-' \
            r'(?P<mode>[HWDPC])' \
            r'(?P<level>1\.[015])' \
            r'(?P<proc>G|_)' \
            r'(?P<proj>[UPML_])' \
            r'(?P<orbit_dir>A|D)$'

ceos_psr2 = r'^LED-ALOS2' \
            r'(?P<orbit>[0-9]{5})' \
            r'(?P<frame>[0-9]{4})-' \
            r'(?P<date>[0-9]{6})-' \
            r'(?P<mode>SBS|UBS|UBD|HBS|HBD|HBQ|FBS|FBD|FBQ|WBS|WBD|WWS|WWD|VBS|VBD)' \
            r'(?P<look_dir>L|R)' \
            r'(?P<level>1\.0|1\.1|1\.5|2\.1|3\.1)' \
            r'(?P<proc>[GR_])' \
            r'(?P<proj>[UPML_])' \
            r'(?P<orbit_dir>A|D)$'

eorc_psr = r'^PSR2-' \
           r'(?P<prodlevel>SLTR)_' \
           r'(?P<pathnr>RSP[0-9]{3})_' \
           r'(?P<date>[0-9]{8})' \
           r'(?P<mode>FBD|WBD)' \
           r'(?P<beam>[0-9]{2})' \
           r'(?P<orbit_dir>A|D)' \
           r'(?P<look_dir>L|R)_' \
           r'(?P<replay_id1>[0-9A-Z]{16})-' \
           r'(?P<replay_id2>[0-9A-Z]{5})_' \
           r'(?P<internal>[0-9]{3})_' \
           r'HDR$'

esa = r'(?P<product_id>(?:SAR|ASA)_(?:IM(?:S|P|G|M|_)|AP(?:S|P|G|M|_)|WV(?:I|S|W|_)|WS(?:M|S|_))_[012B][CP])' \
      r'(?P<processing_stage_flag>[A-Z])' \
      r'(?P<originator_ID>[A-Z\-]{3})' \
      r'(?P<start_day>[0-9]{8})_' \
      r'(?P<start_time>[0-9]{6})_' \
      r'(?P<duration>[0-9]{8})' \
      r'(?P<phase>[0-9A-Z]{1})' \
      r'(?P<cycle>[0-9]{3})_' \
      r'(?P<relative_orbit>[0-9]{5})_' \
      r'(?P<absolute_orbit>[0-9]{5})_' \
      r'(?P<counter>[0-9]{4,})\.' \
      r'(?P<satellite_ID>[EN][12])'

safe = r'^(?P<sensor>S1[ABCD])_' \
       r'(?P<beam>S1|S2|S3|S4|S5|S6|IW|EW|WV|EN|N1|N2|N3|N4|N5|N6|IM)_' \
       r'(?P<product>SLC|GRD|OCN)' \
       r'(?P<resolution>F|H|M|_)_' \
       r'(?P<processingLevel>1|2)' \
       r'(?P<category>S|A)' \
       r'(?P<pols>SH|SV|DH|DV|VV|HH|HV|VH)_' \
       r'(?P<start>[0-9]{8}T[0-9]{6})_' \
       r'(?P<stop>[0-9]{8}T[0-9]{6})_' \
       r'(?P<orbitNumber>[0-9]{6})_' \
       r'(?P<dataTakeID>[0-9A-F]{6})_' \
       r'(?P<productIdentifier>[0-9A-F]{4})' \
       r'\.SAFE$'

tsx = r'^(?P<sat>T[DS]X1)_SAR__' \
      r'(?P<prod>SSC|MGD|GEC|EEC)_' \
      r'(?P<var>____|SE__|RE__|MON1|MON2|BTX1|BRX2)_' \
      r'(?P<mode>SM|SL|HS|HS300|ST|SC)_' \
      r'(?P<pols>[SDTQ])_' \
      r'(?:SRA|DRA)_' \
      r'(?P<start>[0-9]{8}T[0-9]{6})_' \
      r'(?P<stop>[0-9]{8}T[0-9]{6})(?:\.xml|)$'

tdm = r'^(?P<sat>T[D]M1)_SAR__' \
      r'(?P<prod>COS)_' \
      r'(?P<var>____|MONO|BIST|ALT1|ALT2)_' \
      r'(?P<mode>SM|SL|HS)_' \
      r'(?P<pols>[SDQ])_' \
      r'(?:SRA|DRA)_' \
      r'(?P<start>[0-9]{8}T[0-9]{6})_' \
      r'(?P<stop>[0-9]{8}T[0-9]{6})(?:\.xml|)$'
