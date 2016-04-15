import psycopg2, os
from glob import glob

try:
    conn = psycopg2.connect("dbname='GEO413' user='user' host='localhost' password='user'")
except:
    print "database connection failed"

cur = conn.cursor()

try:
    cur.execute("""CREATE TABLE temp(
        Stations_ID int,
        Date date,
        Quality int,
        AIRTEMP float,
        VAPORPRESSURE float,
        COVERAGE float,
        AIRPRESSURE float,
        REL_HUM float,
        V_WIND float,
        AIRTEMP_MAX float,
        AIRTEMP_MIN float,
        AIRTEMP_GROUND_MIN float,
        WIND_MAX float,
        PREC float,
        PREC_IND int,
        SUNHOUR float,
        SNOW int,
        eor varchar(3)
        )""")
except:
    print "no table created"
    conn.rollback()
conn.commit()

# upload DWD textfiles into table 'temp'
path = "/home/user/Desktop/exchange/DWD"
for root, dirs, files in os.walk(path):
    for name in files:
        if name.startswith("produkt") and name.endswith(".txt"):
            product = open(os.path.join(root, name), 'r')
            try:
                cur.copy_expert(sql="COPY temp FROM STDIN with csv header delimiter ';'", file=product)
            except:
                print "upload error for file", name
                conn.rollback()
            conn.commit()
            product.close()
cur.close()
conn.close()
