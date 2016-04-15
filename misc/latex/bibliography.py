
from os import rename, remove

bblname = "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/MasterThesis_Main.bbl"
bibname = "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/Masterarbeit.bib"
txtname = "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/Masterarbeit.txt"

remove(bblname)
remove(bibname)
rename(txtname, bibname)
