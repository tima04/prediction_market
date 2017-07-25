#library(data.table)
## change salanie data format to gandhi format           
    
def main(out_fl = "../../data/salanie_long.csv"):
    data =  open("../../data/salanie_wide.csv").readlines()
    #data = open("/home/kothiyal/Documents/data/data4mJournals/horse_race/salanie/cours_all").readlines()
    columns = ['raceid', 'year', 'winner', 'numhorses', 'winodds']
    sep = "     " 
    raceid = 0
    rslt = ""
    for line in data:
        raceid += 1
        xs = [x.strip() for x in line.strip().split(sep)]
        long = wide2long(raceid, xs)
        if long:
            rslt += long 

    ##import pdb; pdb.set_trace()
    out = open(out_fl, "w")
    out.writelines(",".join(columns) + "\n")
    out.writelines(rslt)
    return None 


def wide2long(raceid, xs):
    rslt = ""
    try:
        year = int(float(xs[0]))
        nhorse = int(float(xs[1]))
        winner = [1] + [0]*(nhorse - 1)
        winodds = xs[2 : (nhorse+2)]
    except:
        return rslt
    for i in range(nhorse):
        rslt += "%s, %s, %s, %s, %s"%(raceid, year, winner[i], nhorse, winodds[i]) + "\n"

    return rslt
        

if __name__ == '__main__':
    main()


