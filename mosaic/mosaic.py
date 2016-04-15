
import raster
import argparse
from imad import imad
import subprocess as sp
from mosaic_aux import *
from time import asctime
from radnorm import radnorm


# helper function for deleting intermediate files
def rollback(files):
    for item in [a+b for a in files for b in ["", ".hdr", ".aux.xml"] if os.path.isfile(a+b)]:
        os.remove(item)


def mosaic(master, slaves, mosaik):
    print "#################################"
    dir_out = os.path.dirname(mosaik)
    ref = raster.Raster(master)

    # #######################################################
    # create a text file containing names of files that have already been included
    includes = mosaik+"_included"
    if not os.path.isfile(includes):
        with open(includes, "a") as infile:
            infile.write(os.path.basename(master)+"\n")
        included = [os.path.basename(master)]
    else:
        with open(includes, "r") as infile:
            included = [line.rstrip() for line in infile]
        print "{0} file{1} been included earlier".format(len(included), " has" if len(included) == 1 else "s have")

    excludes = mosaik+"_excluded"
    if not os.path.isfile(excludes):
        excluded = []
    else:
        with open(excludes, "r") as infile:
            excluded = [line.rstrip() for line in infile]
        print "{0} file{1} been excluded earlier".format(len(excluded), " has" if len(included) == 1 else "s have")

    # reduce slave list to files which have not been included or excluded before
    slaves = [slave for slave in slaves if os.path.basename(slave) not in included and os.path.basename(os.path.splitext(slave)[0])+"_reproject" not in included]
    slaves = [slave for slave in slaves if os.path.basename(slave) not in excluded and os.path.basename(os.path.splitext(slave)[0])+"_reproject" not in excluded]

    # #######################################################
    # gather names of files, which have to be reprojected
    reprojection = {}
    for slave in list(slaves):
        if raster.Raster(slave).proj4 != ref.proj4:
            outname = os.path.join(dir_out, os.path.basename(os.path.splitext(slave)[0])+"_reproject")
            if not os.path.isfile(outname):
                reprojection[slave] = outname
            del slaves[slaves.index(slave)]
            slaves.append(outname)

    # reproject files
    if len(reprojection) > 0:
        print "reprojecting a total of {0} file{1}:".format(len(reprojection), "" if len(reprojection) == 1 else "s"), asctime()
        for item in reprojection:
            print os.path.basename(item)
            raster.reproject(item, ref, reprojection[item])
        print "--------"

    # #######################################################
    # start mosaicing
    while len(slaves) > 0:

        ref = raster.Raster(master)

        # compute intersection areas between master and slaves
        inter_extent = [raster.intersect(ref, raster.Raster(x)) for x in slaves]
        inter_area = [x.area if x is not None else 0 for x in inter_extent]

        if max(inter_area) == 0:
            with open(excludes, "a") as infile:
                for slave in slaves:
                    infile.write(os.path.basename(slave)+"\n")
            return "no intersecting images remaining\n#################################"

        # for all intersecting slaves compute the intersection area of valid pixels
        inter_valid = []
        for i in range(len(slaves)):
            if inter_area[i] > 0:
                rast = raster.Raster(slaves[i])
                mat1 = ref.matrix(1, ref.crop(inter_extent[i]))
                mat2 = rast.matrix(1, rast.crop(inter_extent[i]))
                inter_valid.append(len(mat1[(mat1 != ref.nodata) & (mat2 != rast.nodata)]))
                del mat1, mat2
            else:
                inter_valid.append(0)

        # choose slave based on largest valid pixel overlap with master scene
        index = inter_valid.index(max(inter_valid))
        slave = slaves[index]

        # #######################################################
        print "integrating image {0} ({1} included, {2} excluded, {3} remaining)".format(os.path.basename(slave), len(included), len(excluded), len(slaves))

        slave_base = os.path.basename(os.path.splitext(slave)[0])
        master_base = os.path.basename(os.path.splitext(master)[0])
        slave_sub = os.path.join(dir_out, slave_base+"_sub")
        master_sub = os.path.join(dir_out, master_base+"_sub")
        imad_file = os.path.join(dir_out, master_base+"_"+slave_base+"_imad")
        slave_norm = os.path.join(dir_out, slave_base+"_norm")
        mosaik_temp = mosaik+"_temp"
        deletes = [slave_sub, master_sub, imad_file, slave_norm]
        # #######################################################
        print "--------\n...cropping images to common extent:", asctime()
        rast = raster.Raster(slave)
        rast.crop(ref, slave_sub)
        ref.crop(rast, master_sub)
        # #######################################################
        print "--------\n...detecting changes:", asctime()
        try:
            imad(master_sub, slave_sub, outfile=imad_file)
        except ValueError:
            print "...failed"
            del slaves[index]
            excludes.append(os.path.basename(slave))
            with open(excludes, "a") as infile:
                infile.write(os.path.basename(slave)+"\n")
            rollback(deletes)
            continue
        # #######################################################
        print "--------\n...normalizing images:", asctime()
        try:
            radnorm(master_sub, slave_sub, imad_file, slave, slave_norm)
        except ValueError:
            print "...failed"
            del slaves[index]
            excludes.append(os.path.basename(slave))
            with open(excludes, "a") as infile:
                infile.write(os.path.basename(slave)+"\n")
            rollback(deletes)
            continue
        # #######################################################
        print "--------\n...mosaicing files:", asctime()
        # the order of first slave then master is important as they are stacked on top of each other in the order placed!!!
        sp.check_call(["gdalwarp", "-q", "-overwrite", "-of", "ENVI",
                       "-srcnodata", str(ref.nodata), "-dstnodata", str(ref.nodata),
                       slave_norm, master, mosaik_temp])
        # #######################################################
        # cleaning up
        with open(includes, "a") as infile:
            infile.write(os.path.basename(slave)+"\n")
        included.append(os.path.basename(slave))

        if os.path.isfile(mosaik):
            deletes.append(mosaik)
        if "reproject" in slave:
            deletes.append(slave)

        rollback(deletes)

        os.rename(mosaik_temp, mosaik)
        os.rename(mosaik_temp+".hdr", mosaik+".hdr")
        os.rename(mosaik_temp+".aux.xml", mosaik+".aux.xml")
        raster.Raster(mosaik).png([5, 4, 3], mosaik+"_quicklook.png")
        master = mosaik

        del slaves[index]
        print "#################################"
    return "no more images to include\n#################################"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("master", help="reference image", type=str)
    parser.add_argument("slaves", nargs="+", help="slave image(s) to be mosaiced to the master", type=str)
    parser.add_argument("-o", "outfile", help="name of the mosaic to be created", type=str)
    args = parser.parse_args()
    mosaic(args.master, args.slaves, args.outfile)
