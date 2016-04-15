

import subprocess as sp


def dissolve(inlist):
    out = []
    for i in inlist:
        i = list(i) if type(i) is tuple else i
        out.extend(dissolve(i)) if type(i) is list else out.append(i)
    return out

files = ["D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/MasterThesis_Main.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_conclusion.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_introduction.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_discussion.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_interferometry.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_modeling.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_modeling_AGB.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_modeling_setup.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_modeling_theory.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_modeling_THF.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_orbitmasking.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_polarimetry.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_sampling.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_sar.tex",
         "D:/Dokumente/Studium/Master/GEO512-Masterarbeit/Latex/modules/chapters/MasterThesis_sarprocessing.tex"]

sp.check_call(dissolve(["perl", "D:/TeXcount_3_0_0_24/texcount.pl", "-merge", files]))

# end 1307: 12189/139/970  -> 13298
# end 1407: 13424/139/1063 -> 14626 (+1328)
# end 1507: 13954/170/1094 -> 15218 (+592)
# end 1607: 14906/171/1094 -> 16171 (+953)
# end 1707: 15821/170/1094 -> 17085 (+914)
# end 1807: 16433/178/1201 -> 17812 (+727)
# end 2007: 20297/180/1206 -> 21683 (+3871) ;-)