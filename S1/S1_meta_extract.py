#! /usr/bin/python
# 2015-03-09, Christian Siegert
# modified 2015-08-10, John Truckenbrodt
# This script extracts metadata from Sentinel zip files and writes it to a
# CSV file.
import csv
from glob import glob
import os.path
import re
import xml.etree.ElementTree as ElementTree
import zipfile
import sys


# Regular expressions that match the KML filename and metadata XML filenames in
# a Sentinel zip file.
regExpKmlFile = re.compile("[^/]+/preview/map-overlay\.kml")
regExpMetadataFile = re.compile("[^/]+/annotation/[^/]+\.xml")


# inputDir is the directory where the zip files are located;
# outputFilename is where the CSV file should be saved
def main(inputDir, outputFilename):
    # skipExisting, if true, makes the script skip files that already have an
    # entry in the CSV file.
    skipExisting = True

    # fieldsToExport defines the structure of the CSV file and specifies where
    # the values can be found. It is a list of dictionaries. Each dictionary has
    # a mandatory attribute "csvColumnName" and an optional attribute "xpath".
    # If "xpath" is not set, the CSV column will be created but no value
    # will be written to it. If "xpath" is set, it must be an XPath to the XML
    # element in the XML metadata file, see
    # https://docs.python.org/2/library/xml.etree.elementtree.html#xpath-support
    # Coordinates ("upperLeftLat", "upperLeftLon", ...) are read from the
    # zip file's KML file.
    fieldsToExport = [
        {"csvColumnName": "missionId", "xpath": "./adsHeader/missionId"},
        {"csvColumnName": "productType", "xpath": "./adsHeader/productType"},
        {"csvColumnName": "polarisation", "xpath": "./adsHeader/polarisation"},
        {"csvColumnName": "mode", "xpath": "./adsHeader/mode"},
        {"csvColumnName": "swath", "xpath": "./adsHeader/swath"},
        {"csvColumnName": "startTime", "xpath": "./adsHeader/startTime"},
        {"csvColumnName": "stopTime", "xpath": "./adsHeader/stopTime"},
        {"csvColumnName": "absoluteOrbitNumber", "xpath": "./adsHeader/absoluteOrbitNumber"},
        {"csvColumnName": "missionDataTakeId", "xpath": "./adsHeader/missionDataTakeId"},
        {"csvColumnName": "imageNumber", "xpath": "./adsHeader/imageNumber"},
        {"csvColumnName": "pass", "xpath": "./generalAnnotation/productInformation/pass"},
        {"csvColumnName": "outputPixels", "xpath": "./imageAnnotation/imageInformation/outputPixels"},
        {"csvColumnName": "rangePixelSpacing", "xpath": "./imageAnnotation/imageInformation/rangePixelSpacing"},
        {"csvColumnName": "azimuthPixelSpacing", "xpath": "./imageAnnotation/imageInformation/azimuthPixelSpacing"},
        {"csvColumnName": "numberOfSamples", "xpath": "./imageAnnotation/imageInformation/numberOfSamples"},
        {"csvColumnName": "numberOfLines", "xpath": "./imageAnnotation/imageInformation/numberOfLines"},
        {"csvColumnName": "upperLeftLat"},  # Value is read from KML file
        {"csvColumnName": "upperLeftLon"},  # Value is read from KML file
        {"csvColumnName": "upperRightLat"},  # Value is read from KML file
        {"csvColumnName": "upperRightLon"},  # Value is read from KML file
        {"csvColumnName": "lowerRightLat"},  # Value is read from KML file
        {"csvColumnName": "lowerRightLon"},  # Value is read from KML file
        {"csvColumnName": "lowerLeftLat"},  # Value is read from KML file
        {"csvColumnName": "lowerLeftLon"},  # Value is read from KML file
    ]

    # existingFilenames is a dictionary of existing files in the CSV file. Used
    # as look-up table.
    existingFilenames = {}

    # existingRowCount is the number of rows in the CSV file (inluding header).
    existingRowCount = 0

    # If files with an existing entry in the CSV file should be skipped, compile
    # a filename look-up table.
    if skipExisting:
        filenames, rowCount, err = getFilenamesAndRowCountFromCsv(outputFilename)
        if err != "":
            print "ERROR: %s" % err
            return
        existingFilenames = filenames
        existingRowCount = rowCount

    filenamePattern = os.path.join(inputDir, "*.zip")
    filenames = glob(filenamePattern)
    print "Found %d zip files in %s" % (len(filenames), inputDir)
    print "Extracting metadata..."

    metadataLists = []
    newCount = 0
    skipCount = 0
    errorCount = 0

    for filename in filenames:
        if skipExisting and existsInCsv(existingFilenames, filename):
            skipCount += 1
            continue

        if not zipfile.is_zipfile(filename):
            print "Skipping invalid zip file %s" % filename
            errorCount += 1
            continue

        metadataList, err = getMetadata(filename, fieldsToExport)
        if err != "":
            print "ERROR: %s: %s" % (os.path.basename(filename), err)
            errorCount += 1
            continue

        if len(metadataList) == 0:
            print "ERROR: No metadata file found in %s" % os.path.basename(filename)
            errorCount += 1
            continue

        metadataLists.append(metadataList)
        newCount += 1

    writeCsvFile(outputFilename, skipExisting, existingRowCount + 1, fieldsToExport, metadataLists)
    print "Extracted metadata from %d zip files, skipped %s due to existing CSV entry, skipped %d due to errors." % (newCount, skipCount, errorCount)


def getFilenamesAndRowCountFromCsv(csvFilename):
    filenames = {}
    rowCount = 0

    if not os.path.isfile(csvFilename):
        return filenames, rowCount, ""

    csvFile = open(csvFilename, "r")
    reader = csv.DictReader(csvFile)

    for row in reader:
        if "zipFilename" not in row:
            return None, 0, "CSV file is missing column zipFilename"
        filenames[row["zipFilename"]] = True
        rowCount += 1
    csvFile.close()
    return filenames, rowCount, ""


# existsInCsv returns true if zipFilename exists in the CSV file.
def existsInCsv(existingFilenames, filename):
    return os.path.basename(filename) in existingFilenames


# getMetadata opens the zip file identified by filename and returns a list of
# Metadata instances.
def getMetadata(filename, fieldsToExport):
    zipFile = zipfile.ZipFile(filename, "r")
    metadataList = []
    coordinateFields = []

    # Read coordinates from KML file
    for info in zipFile.infolist():
        if not regExpKmlFile.match(info.filename):
            continue

        xml = zipFile.read(info)
        coordinateFields, err = getCoordinates(xml)

        if err != "":
            return [], err

    # Read metadata from XML files
    for info in zipFile.infolist():
        if not regExpMetadataFile.match(info.filename):
            continue

        xml = zipFile.read(info)
        metadataFields = getMetadataFields(xml, fieldsToExport)
        metadataFields.extend(coordinateFields)

        metadata = Metadata(filename, info.filename, metadataFields)
        metadataList.append(metadata)

    zipFile.close()
    return metadataList, ""


# getCoordinates parses the provided XML of the KML file for coordinates and
# returns a list of MetadataField instances. If an error occurred, the returned
# list is None and the error string contains a message.
def getCoordinates(xml):
    rootElement = ElementTree.fromstring(xml)

    elements = rootElement.findall(".//coordinates")

    if elements is None or len(elements) == 0:
        return None, "Field 'coordinates' not found in KML file."
    elif len(elements) >= 2:
        return None, "Found field 'coordinates' multiple times in KML file."

    coordinatePairs = elements[0].text.split(" ")
    if len(coordinatePairs) != 4:
        return None, "expected 4 coordinate pairs, got %d" % len(coordinatePairs)

    coordinates = []

    for coordinatePair in coordinatePairs:
        pieces = coordinatePair.split(",")
        if len(pieces) < 2:
            return None, "expected 2 coordinates in coordinate pair, got %d" % len(pieces)
        coordinates.append((float(pieces[1]), float(pieces[0])))

    latLonQuad = LatLonQuad(coordinates[0], coordinates[1], coordinates[2], coordinates[3])
    latLonQuad.sort()

    coordinateFields = [
        MetadataField("upperLeftLat", latLonQuad.upperLeft[0]),
        MetadataField("upperLeftLon", latLonQuad.upperLeft[1]),
        MetadataField("upperRightLat", latLonQuad.upperRight[0]),
        MetadataField("upperRightLon", latLonQuad.upperRight[1]),
        MetadataField("lowerRightLat", latLonQuad.lowerRight[0]),
        MetadataField("lowerRightLon", latLonQuad.lowerRight[1]),
        MetadataField("lowerLeftLat", latLonQuad.lowerLeft[0]),
        MetadataField("lowerLeftLon", latLonQuad.lowerLeft[1]),
    ]

    return coordinateFields, ""


# parse XML of the metadata file and return a list of MetadataField instances
def getMetadataFields(xml, fieldsToExport):
    metadataFields = []
    rootElement = ElementTree.fromstring(xml)

    for fieldToExport in fieldsToExport:
        if "xpath" not in fieldToExport:
            continue

        elements = rootElement.findall(fieldToExport["xpath"])

        if len(elements) == 0:
            print "WARNING: Field %s is missing in %s, metadata file: %s" % (fieldToExport["csvColumnName"], item["zipFilename"], csvRow["metadataFilename"])
            continue
        elif len(elements) >= 2:
            print "WARNING: Metadata contains more than one element %s in %s, metadata file: %s" % (fieldToExport["csvColumnName"], item["zipFilename"], csvRow["metadataFilename"])
            continue

        metadataField = MetadataField(fieldToExport["csvColumnName"], elements[0].text)
        metadataFields.append(metadataField)
    return metadataFields


# writeCsvFile writes metadataLists to outputFilename in CSV format.
def writeCsvFile(outputFilename, skipExisting, startId, fieldsToExport, metadataLists):
    csvColumnNames = ["id", "zipFilename", "metadataFilename"]

    for fieldToExport in fieldsToExport:
        csvColumnNames.append(fieldToExport["csvColumnName"])

    csvFileExists = os.path.exists(outputFilename)

    fileMode = "w"
    if skipExisting: fileMode = "a"

    outputFile = open(outputFilename, fileMode)
    writer = csv.DictWriter(outputFile, csvColumnNames)

    if not csvFileExists or not skipExisting:
        writer.writeheader()

    id = startId

    for metadataList in metadataLists:
        for metadata in metadataList:
            csvRow = {
                "id": id,
                "zipFilename": os.path.basename(metadata.zipFilename),
                "metadataFilename": os.path.basename(metadata.filename),
            }

            for field in metadata.fields:
                csvRow[field.csvColumnName] = field.value

            writer.writerow(csvRow)
            id += 1
    outputFile.close()


# LatLonQuad represents a quad where each vertex is a geographical coordinate.
class LatLonQuad:
    # Each coordinate is a tuple of the form (<lat>, <lon>).
    def __init__(self, coordinate1, coordinate2, coordinate3, coordinate4):
        self.coordinates = [
            coordinate1,
            coordinate2,
            coordinate3,
            coordinate4,
        ]

        # Default values
        self.lowerLeft = (0, 0)
        self.lowerRight = (0, 0)
        self.upperLeft = (0, 0)
        self.upperRight = (0, 0)

    # sort sorts the coordinates by their position (upper left, upper right,
    # lower right, lower left).
    def sort(self):
        sortedByLat = self.sortByLat(self.coordinates)
        upper = self.sortByLon(sortedByLat[:2])
        lower = self.sortByLon(sortedByLat[2:])

        self.lowerLeft = lower[0]
        self.lowerRight = lower[1]
        self.upperLeft = upper[0]
        self.upperRight = upper[1]

    # sortByLat returns a copy of coordinates sorted by latitude.
    def sortByLat(self, coordinates):
        return sorted(coordinates, key=lambda coordinate: coordinate[0], reverse=True)

    # sortByLon returns a copy of coordinates sorted by longitude.
    def sortByLon(self, coordinates):
        return sorted(coordinates, key=lambda coordinate: coordinate[1])


# Metadata contains the data of a zip file's metadata file.
class Metadata:
    def __init__(self, zipFilename, filename, fields):
        self.fields = fields
        self.filename = filename
        self.zipFilename = zipFilename


# MetadataField contains the value that should be written into the specified
# CSV column.
class MetadataField:
    def __init__(self, csvColumnName, value):
        self.csvColumnName = csvColumnName
        self.value = value


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
