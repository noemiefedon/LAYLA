import os
import fnmatch

def findReplace(directory, find, replace, filePattern):
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            with open(filepath) as f:
                s = f.read()
            s = s.replace(find, replace)
            with open(filepath, "w") as f:
                f.write(s)

find = "from src.LAYLA_V02.panels import Panel"
replace = "from src.LAYLA_V02.panels import Panel"

findReplace(r'C:\LAYLA_and_BELLA', find, replace, '*.py')
findReplace(r'C:\LAYLA_and_BELLA\RELAY', find, replace, '*.py')
