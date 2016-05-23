'''Script to initially fill the database of a folder
'''
from pyroSAR import identify
from ancillary import finder
import getpass


def fill_db(folder):
    scenes = finder(folder, ["^(?:S1[AB]|ASA|SAR)"], regex=True)
    print scenes
    print getpass.getuser()
    selection = []
    print "collecting files"
    for scene in scenes:
        try:
            id = identify(scene)

        except IOError:
            continue
        #print os.path.dirname(id.file)
        #print os.path.dirname(id.scene)
        try:
            id.export2sqlite()
        except RuntimeError as e:
            print e

fill_db('../data/')
