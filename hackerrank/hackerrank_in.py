
import subprocess as sp

with open("E:/hackerrank/testcase.txt", "r") as infile:
    testcases = infile.read()

out, err = sp.Popen(["python", "sherlockbeast.py"], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate(testcases)

print out, err
