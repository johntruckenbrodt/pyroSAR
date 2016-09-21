import os
import zipfile
from ftplib import FTP

pathlist = ['p' + format(x, '03d') for x in range(115, 124)]
rowlist = ['r' + format(x, '03d') for x in range(23, 34)]

path_local = 'D:/Landsat/VCF'
path_server = 'ftp.glcf.umd.edu'
path_remote = 'glcf/LandsatTreecover/WRS2'

ftp = FTP(path_server)

ftp.login()
for path in pathlist:
    for row in rowlist:
        path_target = os.path.join(path_remote, path, row, path + row + '_TC_2000')
        try:
            ftp.cwd(path_target)
            file_target = path + row + '_TC_2000.tif.gz'
            file_local = os.path.join(path_local, file_target)
            # ftp.retrlines('LIST')
            # ftp.nlst()
            if not os.path.isfile(file_local):
                print file_target, '->', file_local
                ftp.retrbinary('RETR ' + file_target, open(file_local, 'wb').write)
                with zipfile.ZipFile(file_local, 'r') as z:
                    z.extractall(path_local)
                os.remove(file_local)
        except:
            continue
ftp.quit()
