import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="myFile", help="Open specified file")
args = parser.parse_args()
myFile = args.myFile

elementsToDelete = 0

f = open(myFile,'r')

for idx, line in enumerate(f.readlines()):
    idx += 1
    if '$EndElements' in line.strip():
        finalElementIdx = idx-1
        break
    if '$Elements' in line.strip():
        maxElementsIdx = idx
        firstElementIdx = idx+1

f.close()

f = open(myFile,'r')
lines = f.readlines()
maxElements = int(lines[maxElementsIdx].strip())

els_list = []
idx_list = 0
new_el_id = 1
for i in range(int(maxElements)):
    temp_list = lines[i+firstElementIdx].strip().split()
    
    if (temp_list[1] == "15" or temp_list[1] == "1" or temp_list[1] == "26"):
        elementsToDelete += 1
    else:
        temp_list[3] = temp_list[4]
        temp_list[0] = str(new_el_id)
        new_el_id += 1
        temp_list = " ".join(temp_list)
        temp_list = temp_list + "\n"
        els_list.append(temp_list)

    
    idx_list += 1

f.close()

fr = open(myFile,"r")
fr_lines = fr.readlines()
for idx, line in enumerate(fr_lines):
    if idx in range(firstElementIdx+elementsToDelete, finalElementIdx):
        fr_lines[idx] = els_list[idx-(firstElementIdx+elementsToDelete)]

fr.close()

fr = open(myFile,"r")
fr_lines = fr.readlines()
lines_to_write = []
for idx, line in enumerate(fr_lines):
    if idx in range(0, firstElementIdx):
        if idx == firstElementIdx-1:
            lines_to_write.append(str(maxElements-elementsToDelete)+"\n")
        else:
            lines_to_write.append(fr_lines[idx])
    if idx in range(firstElementIdx+elementsToDelete, finalElementIdx):
        lines_to_write.append(els_list[idx-firstElementIdx-elementsToDelete])

lines_to_write.append("$EndElements"+"\n")

fr.close()

# os.remove(myFile)

fw = open(myFile,"w")
fw.writelines(lines_to_write)
fw.close()

print("Mesh converted to Adapter compatible.")
                