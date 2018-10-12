lines = open(readlines, "FileCG.txt")
newlines = sort(lines)
writedlm("FileCG.txt",newlines)

lines = open(readlines, "FileCR.txt")
newlines = sort(lines)
writedlm("FileCR.txt",newlines)
