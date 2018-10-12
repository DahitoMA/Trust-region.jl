fichier = open("NewResults_LS/tabLSQR.txt", "r")
file = open("NewResults_LS/tabLSQR_bis.txt", "w")
for line in eachline(fichier)
	elements = split(line," ")
	for i=1:length(elements) - 4
    	write(file, elements[i])
    	write(file, " ")
    end
    write(file, "\n")
end
close(file)
close(fichier)
