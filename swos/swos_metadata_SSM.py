import os
import re
import xml.etree.ElementTree as ET
from calendar import monthrange
from time import ctime, strftime, strptime

from ancillary import finder
from spatial import raster
from xml_util import getNamespaces

template = "/geonfs01_vol1/ve39vem/swos_test/meta_xml/template_SSM.xml"

namespaces = getNamespaces(template)

sites = ["Egypt_Burullus", "France_Camargue", "Greece_Eastern-Macedonia", "Jordan_Azraq", "Spain_Fuente-de-Piedra", "Tanzania_Kilombero", "Sweden_Store-Mosse"]

update = True

for site in sites:
    print site
    files = finder("/geonfs01_vol1/ve39vem/swos_process/{}".format(site), ["SSM_[0-9_]+\.tif"], regex=True)
    for item in files:
        outname = item.replace(".tif", ".xml")

        if not update and os.path.isfile(outname):
            continue

        with open(template, "r") as infile:
            tree = ET.fromstring(infile.read())

        name = tree.find('.//ns0:fileIdentifier/./', namespaces)

        ###########################################
        meta_contact = tree.find('.//ns0:contact', namespaces)
        meta_contact_organization = meta_contact.find('.//ns0:organisationName/./', namespaces)
        meta_contact_email = meta_contact.find('.//ns0:electronicMailAddress/./', namespaces)
        ###########################################
        meta_date = tree.find('.//ns0:dateStamp/ns2:Date', namespaces)
        ###########################################
        meta_srs = tree.find('.//ns0:referenceSystemInfo/', namespaces)
        meta_srs_url = meta_srs.find('.//ns0:RS_Identifier/ns0:code/./', namespaces)
        ###########################################
        data = tree.find('.//ns0:MD_DataIdentification', namespaces)
        data_title = data.find('.//ns0:title/./', namespaces)

        data_date = data.find('.//ns0:date/ns2:Date', namespaces)
        data_date_type = data.find('.//ns0:dateType/ns0:CI_DateTypeCode', namespaces)
        data_date_value = data.find('.//ns0:dateType/ns0:CI_DateTypeCode', namespaces).text

        data_abstract = data.find('.//ns0:abstract/./', namespaces)

        data_contact_organization = data.find('.//ns0:pointOfContact/.//ns0:organisationName/./', namespaces)
        data_contact_email = data.find('.//ns0:electronicMailAddress/./', namespaces)
        ####################
        keywords = data.findall('.//ns0:descriptiveKeywords', namespaces)
        # for keyword in keywords:
        #     print keyword.find('.//ns0:MD_Keywords/ns0:keyword/./', namespaces)
        ####################

        data_extent = data.find('.//ns0:extent/ns0:EX_Extent', namespaces)

        data_extent_bbox = data_extent.find('.//ns0:geographicElement/ns0:EX_GeographicBoundingBox', namespaces)
        data_extent_bbox_west = data_extent_bbox.find('.//ns0:westBoundLongitude/./', namespaces)
        data_extent_bbox_east = data_extent_bbox.find('.//ns0:eastBoundLongitude/./', namespaces)
        data_extent_bbox_north = data_extent_bbox.find('.//ns0:northBoundLatitude/./', namespaces)
        data_extent_bbox_south = data_extent_bbox.find('.//ns0:southBoundLatitude/./', namespaces)

        data_res = data.find('.//ns0:spatialResolution/.//ns2:Distance', namespaces)

        data_extent_date_begin = data_extent.find('.//ns0:temporalElement/ns0:EX_TemporalExtent/ns0:extent/ns3:TimePeriod/ns3:beginPosition', namespaces)
        data_extent_date_end = data_extent.find('.//ns0:temporalElement/ns0:EX_TemporalExtent/ns0:extent/ns3:TimePeriod/ns3:endPosition', namespaces)

        #######################################################################################################################################################

        name.text = os.path.basename(os.path.splitext(item)[0])

        data_title.text = "SWOS Surface Soil Moisture"

        meta_contact_organization.text = "Friedrich-Schiller-University Jena"
        meta_contact_email.text = "john.truckenbrodt@uni-jena.de"

        meta_date.text = strftime("%Y-%m-%d")

        data_abstract.text = "This file contains an monthly means of a SAR-derived relative surface soil moisture index. " \
                             "Values are scaled in percent from 0 to 100. Data gaps are marked with value 255. " \
                             "The data set was created in the framework of SWOS (Satellite-based Wetland Observation System)."

        data_contact_organization.text = "Friedrich-Schiller-University Jena"
        data_contact_email.text = "john.truckenbrodt@uni-jena.de"

        ras = raster.Raster(item)
        bbox = ras.bbox()
        bbox.reproject(meta_srs_url.text)
        ext = bbox.extent
        data_extent_bbox_west.text = str(ext["xmin"])
        data_extent_bbox_east.text = str(ext["xmax"])
        data_extent_bbox_north.text = str(ext["ymax"])
        data_extent_bbox_south.text = str(ext["ymin"])

        data_res.text = str(ras.res[0])
        data_res.attrib['uom'] = "meter"

        try:
            data_date.text = strftime("%Y-%m-%d", strptime(ras.raster.GetMetadataItem("TIFFTAG_DATETIME"), "%Y:%m:%d %H:%M:%S"))
        except TypeError:
            data_date.text = strftime("%Y-%m-%d", strptime(ctime(os.path.getmtime(item)), "%a %b %d %H:%M:%S %Y"))

        data_date_type.text = "creation"
        data_date_type.attrib['codeListValue'] = "creation"

        date = re.search("[0-9_]{7}\.tif", os.path.basename(item)).group()
        year, month = date.strip(".tif").split("_")
        data_extent_date_begin.text = year+"-"+month+"-01"
        data_extent_date_end.text = year+"-"+month+"-{}".format(monthrange(int(year), int(month))[1])

        keywords[1].find(".//ns0:MD_Keywords/ns0:keyword/ns2:CharacterString", namespaces).text = "surface soil moisture"

        with open(outname, "w") as outfile:
            outfile.write(ET.tostring(tree))
