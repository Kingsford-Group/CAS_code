import sys

def stripN ( filename ):
    read_file = open(filename, "r")
    write_file = open(filename + "_noN", "w")
    
    for line in read_file.xreadlines():
        seed = line.split(' ')[0]

        if 'N' not in seed:
            write_file.write(line)

    read_file.close()
    write_file.close()


for arg in sys.argv[1:]:
    print("processing " + arg + " file")
    stripN(arg)
    print("done")
