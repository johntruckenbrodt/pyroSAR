import psycopg2, os
from glob import glob

try:
    conn = psycopg2.connect("dbname='GEO413' user='user' host='localhost' password='user'")
except:
    print "database connection failed"
cur = conn.cursor()
try:
    cur.execute("""CREATE TABLE data_tlug(
        Stations_ID varchar(7),
        date timestamp,
        value float,
        parameter varchar(15),
        quality varchar(10)
        )""")
except:
    print "no table created"
    conn.rollback()
cur.execute("""SET datestyle = "ISO, DMY";""")
conn.commit()

path = "/media/user/Data/Dokumente/Studium/Master/GEO413-Geodatenbanken/GEO413_data/TLUG"
for root, dirs, files in os.walk(path):
    for name in files:
        if name.endswith("_28.txt"):
            product = open(os.path.join(root, name), 'r')
            lines = product.readlines()
            parameter = lines[0].split(" ")[0]
            ID = lines[3][16:].strip('\r\n')
            quality = lines[4].split(",")[0][10:]
            temp = open("temp.txt", "wr")
            lines = lines[5:]
            replace1 = {" keine Werte": "-999"}
            replace2 = {".": "/"}
            for line in lines:
                for key in replace1.keys():
                    line = line.replace(key, replace1[key])
                for key in replace2.keys():
                    rep = line.replace(key, replace2[key])[0:10]
                line.strip('\r\n')
                line = ID + ", " + line.strip('\r\n') + ", " + parameter + ", " + quality + '\r\n'
                temp.write(line)
            temp.close()
            product.close()
            cur.copy_expert(sql="COPY data_tlug FROM STDIN with delimiter ','", file=open("temp.txt", 'r'))
            conn.commit()
            os.remove("temp.txt")
cur.close()
conn.close()
