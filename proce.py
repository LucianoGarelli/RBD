def proce(file_in):

    num_lines = sum(1 for line in open(file_in, 'r'))
    ff = open(file_in, 'r')
    f_name=file_in.split('.txt')
    f_save = f_name[0]+str('_proc.txt')
    # Lectura comentario
    line = ff.readline()
    line1 = ff.readline()
    frags = line1.split(',')
    V1 = float(frags[0])
    ff2 = open(f_save, 'w')
    ff2.write(line)
    for k in range(num_lines-1):
        line2 = ff.readline()
        frags2 = line2.split(',')
        if not line2:
            ff.seek(0)
            last_line = ff.readlines()[-1]
            ff2.write(last_line)
            break
        V2 = float(frags2[0])
        if V1 != V2 or ff is None:
            ff2.write(line1)
            V1 = V2
        line1 = line2
    ff.close()
    ff2.close()
    return(1)
