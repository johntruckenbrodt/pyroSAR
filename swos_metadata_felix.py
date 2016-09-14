# -*- coding: utf-8 -*-
"""
Metadata Handling for SWOS.
@author Felix Cremer
felix.cremer@uni-jena.de
Should not be used by itself.
"""
import logging
import xml.etree.ElementTree as ET
from ancillary import finder
import os.path
import requests
from vector import Vector


def add_linkage(xml, ext, url='http://artemis.geogr.uni-jena.de/geoserver/', jpeg_folder='http://artemis.geogr.uni-jena.de/ec/swos/products/'):
    """
    Adds the links of the geoserver at url into the xml metadata file.
    :param xml: path to the metadata file
    :param ext: Extension of the underlying product, defines if WCS or WFS is used
    :param url: url to the geoserver
    :return: None
    """
    # Open the xml file
    try:
        tree = ET.parse(xml)
    except IOError as e:
        print e
        return e
    # Set the type information for the OGC download service, from the file extension
    if ext == '.shp':
        type = 'wfs'
        type_name = 'Feature'
        type_spec = '1.1.0'
    elif ext in ['.tif', '.tiff', '.geotiff', '.geotif']:
        type = 'wcs'
        type_name = 'Coverage'
        type_spec = '1.0.0'
    else:
        print 'The File type is not a Raster or Vector type'

    root = tree.getroot()
    # Find the onLine node of the xml in which the links will be integrated.
    online = root.find('.//{http://www.isotc211.org/2005/gmd}onLine')
    # Get the layername from the metadata
    name = root.find('.//{http://www.isotc211.org/2005/gmd}fileIdentifier/{http://www.isotc211.org/2005/gco}CharacterString').text
    # Get the name of the testsite. Works only with the right folder structure
    testsite = os.path.dirname(xml).split('/')[-2]
    # Prepare the information, which will be included into the xml templates.
    testsite_url = ''.join([url, testsite, '/'])
    wms_url = ''.join([testsite_url, 'wms?'])
    type_url = ''.join([testsite_url, type, '?'])
    jpeg_link = '{folder}{jpeg}.jpg'.format(folder=jpeg_folder, jpeg=name)
    data_id = root.find('.//{http://www.isotc211.org/2005/gmd}MD_DataIdentification')
    png = 'png/{z}/{x}/{y}.png'
    tms_url = '{url}gwc/service/tms/1.0.0/{testsite}:{layername}@EPSG:900913@{png}'
    tms_url = tms_url.format(testsite=testsite, layername=name, png=png, url=url)

    # Define the xml templates for the different linkage information as strings
    jpeg_string = '''
    <gmd:graphicOverview xmlns:gmd="http://www.isotc211.org/2005/gmd" xmlns:gco="http://www.isotc211.org/2005/gco">
        <gmd:MD_BrowseGraphic>
          <gmd:fileName>
            <gco:CharacterString>{link}</gco:CharacterString>
          </gmd:fileName>
          <gmd:fileDescription>
            <gco:CharacterString>Thumbnail of the {name} map.</gco:CharacterString>
          </gmd:fileDescription>
          <gmd:fileType>
            <gco:CharacterString>JPEG</gco:CharacterString>
          </gmd:fileType>
        </gmd:MD_BrowseGraphic>
      </gmd:graphicOverview>
    '''

    wms_string = '''
    <gmd:CI_OnlineResource xmlns:gmd="http://www.isotc211.org/2005/gmd" xmlns:gco="http://www.isotc211.org/2005/gco">
            <gmd:linkage>
                <gmd:URL>{url}</gmd:URL>
            </gmd:linkage>
            <gmd:protocol>
                <gco:CharacterString>urn:ogc:serviceType:WebMapService:1.1.1:HTTP</gco:CharacterString>
            </gmd:protocol>
            <gmd:name>
                <gco:CharacterString>{layername}</gco:CharacterString>
            </gmd:name>
            <gmd:description>
                <gco:CharacterString>WMS of the Satellite based Wetlands Observation Services(SWOS) project.</gco:CharacterString>
            </gmd:description>
            <gmd:function>
                <gmd:CI_OnLineFunctionCode
                    codeList="http://www.isotc211.org/2005/resources/Codelist/gmxCodelists.xml#CI_OnLineFunctionCode" codeListValue="download">download</gmd:CI_OnLineFunctionCode>
            </gmd:function>
        </gmd:CI_OnlineResource>
    '''
    tms_string = '''
    <gmd:CI_OnlineResource xmlns:gmd="http://www.isotc211.org/2005/gmd" xmlns:gco="http://www.isotc211.org/2005/gco">
        <gmd:linkage>
            <gmd:URL>{url}</gmd:URL>
        </gmd:linkage>
        <gmd:protocol>
            <gco:CharacterString>urn:essi:serviceType:TiledMapService:1.0.0:HTTP</gco:CharacterString>
        </gmd:protocol>
        <gmd:name>
            <gco:CharacterString>{layername}</gco:CharacterString>
        </gmd:name>
        <gmd:description>
            <gco:CharacterString>TMS of the Satellite based Wetlands Observation Services(SWOS) project.</gco:CharacterString>
        </gmd:description>
        <gmd:function>
            <gmd:CI_OnLineFunctionCode
                codeList="http://www.isotc211.org/2005/resources/Codelist/gmxCodelists.xml#CI_OnLineFunctionCode" codeListValue="browseGraphic">browseGraphic</gmd:CI_OnLineFunctionCode>
        </gmd:function>
    </gmd:CI_OnlineResource>
    '''
    down_string = '''
    <gmd:CI_OnlineResource xmlns:gmd="http://www.isotc211.org/2005/gmd" xmlns:gco="http://www.isotc211.org/2005/gco">
            <gmd:linkage>
                <gmd:URL>{url}</gmd:URL>
            </gmd:linkage>
            <gmd:protocol>
                <gco:CharacterString>urn:ogc:serviceType:Web{type_name}Service:{type_spec}:HTTP</gco:CharacterString>
            </gmd:protocol>
            <gmd:name>
                <gco:CharacterString>{layername}</gco:CharacterString>
            </gmd:name>
            <gmd:description>
                <gco:CharacterString> Web {type_name} Service of the Satellite based Wetlands Observation Services(SWOS) project.</gco:CharacterString>
            </gmd:description>
            <gmd:function>
                <gmd:CI_OnLineFunctionCode
                    codeList="http://www.isotc211.org/2005/resources/Codelist/gmxCodelists.xml#CI_OnLineFunctionCode" codeListValue="download">download</gmd:CI_OnLineFunctionCode>
            </gmd:function>
        </gmd:CI_OnlineResource>
    '''

    # Format the strings for integration
    wms_string = wms_string.format(url=wms_url, layername=name)
    tms_string = tms_string.format(testsite=testsite, layername=name, png=png, url=tms_url)
    down_string = down_string.format(testsite=testsite, layername=name, url=type_url, type_name=type_name, type_spec=type_spec)
    jpeg_string = jpeg_string.format(link=jpeg_link, name=name)

    # Add WMS to xml tree
    wms = ET.XML(wms_string)
    online.append(wms)
    # Add TMS to xml tree
    tms = ET.XML(tms_string)
    online.append(tms)
    # Add Download-Info to xml tree
    download = ET.XML(down_string)
    online.append(download)

    # todo: Add better way to check, if the jpeg is available
    if 'GW2' in name:
        jpeg = ET.XML(jpeg_string)
        data_id.append(jpeg)
    # Write the data to the xml file
    tree.write(xml)
    # Format the data as dictionary, which should be integrated into the Django framework
    data = {'ogc_link': wms_url, 'download_url': type_url}
    log.debug(data)
    return data


# todo: Does it need to be an extra function
def upload_cswt(xml, url='http://localhost:8000/'):
    """
    Upload the Metadata in a xml file into the pycsw instance at url.
    :param xml: The path to the metadata file
    :param url: The pycsw instance
    :return: None
    """
    # Frame of the CSW Transaction
    input = '''
    <csw:Transaction xmlns:csw="http://www.opengis.net/cat/csw/2.0.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.opengis.net/cat/csw/2.0.2 http://schemas.opengis.net/csw/2.0.2/CSW-publication.xsd" service="CSW" version="2.0.2">
    <csw:Insert>
    {metadata}
    </csw:Insert></csw:Transaction>
    '''
    with open(xml) as f:
        metadata = f.read()
        input = input.format(metadata=metadata)
        # todo: Why doesn't it need authentification?
        r = requests.post(url, input)


def delete_duplicate_link(xml):
    """
    Delete the links, which are integrated by add_linkage from the xml file.
    Only used for debugging.
    :param xml: The Metadata file
    :return: None
    """
    tree = ET.parse(xml)
    root = tree.getroot()
    online = root.find('.//{http://www.isotc211.org/2005/gmd}onLine')
    for i, resource in enumerate(online):
        try:
            protocol = resource[1][0].text
        except IndexError:
            protocol = 'link'
        # print protocol
        if protocol != 'link':
            online.remove(resource)
            print 'print yo'
    # print online
    data_id = root.find('.//{http://www.isotc211.org/2005/gmd}MD_DataIdentification')
    # print data_id
    graphic = data_id.find('.//{http://www.isotc211.org/2005/gmd}graphicOverview')
    # print graphic
    try:
        data_id.remove(graphic)
    except ValueError as e:
        pass

    tree.write(xml)


def metadata_from_template(shapefile, template_xml):
    """
    Make a metadata file for shapefile from a Inspire conform template.
    It is assumed, that the shapefile is in the folder-structure /country/wetland/shapefile.shp
    :param shapefile: path to the shapefile
    :param template_xml: path to the template
    :return: None
    """

    shape_name = os.path.splitext(os.path.basename(shapefile))[0]
    # Open the xml file
    template_tree = ET.parse(template_xml)
    template_root = template_tree.getroot()
    # Set the shapefile name as fileIdentifier and MD DataIdentification
    name = template_root.find('.//{http://www.isotc211.org/2005/gmd}fileIdentifier')[0]
    name.text = shape_name
    CI_Citation = template_root.find(".//{http://www.isotc211.org/2005/gmd}identificationInfo/{http://www.isotc211.org/2005/gmd}MD_DataIdentification")[0][0]
    title = CI_Citation[0]
    title[0].text = shape_name
    code = CI_Citation.find('{http://www.isotc211.org/2005/gmd}identifier/{http://www.isotc211.org/2005/gmd}RS_Identifier')[0]
    code[0].text = shape_name
    # Open the shapefile to get the bounding box corners
    shape = Vector(shapefile)
    extent = shape.extent

    # Add the bounding box corners to the xml file
    bounding_box = template_root.find('.//{http://www.isotc211.org/2005/gmd}EX_GeographicBoundingBox')

    bounding_box[0][0].text = str(extent['xmin'])
    bounding_box[1][0].text = str(extent['xmax'])
    bounding_box[2][0].text = str(extent['ymin'])
    bounding_box[3][0].text = str(extent['ymax'])

    # Save the xml tree in shapefile.xml
    xml_path = os.path.splitext(shapefile)[0] + '.xml'
    template_tree.write(xml_path)


if __name__ == '__main__':
    # Delete duplicates, the other functionality should be accessed as imports
    logging.basicConfig(level=logging.DEBUG)
    logging.info('Started')
    folder = '/home/qe89hep/data/SWOS/wetlands'
    metadata_files = finder(folder, ['*.xml'])
    template = '/home/qe89hep/data/SWOS/template_watersheds.xml'
    for metadata in metadata_files:
        logging.debug(metadata)
        # metadata_from_template(metadata, template)
        delete_duplicate_link(metadata)
        # add_linkage(metadata, '.shp')
        # upload_metadata(metadata, '/home/qe89hep/pycsw/default.cfg')
    logging.info('End')
