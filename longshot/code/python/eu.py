import pandas as pd
import numpy as np

#data_fl = "../../data/test_data.csv" 
data_fl = "../../data/dataset_full.csv" 
df = pd.read_csv(data_fl, sep =",")

def cara(xs, theta):
    return (1 - np.exp(-theta*xs))/theta


def solve_ps(t, Rs, epsilon=10e-10):
    n = len(Rs)
    utils = cara(Rs, t)
    u_min_1 = cara(-1, t)
    w = 1/np.sum(1/(utils - u_min_1)) + u_min_1 
    ps = (w - u_min_1)/(utils - u_min_1)
    return ps

def loglike(t, dt):
    #ps = dt.groupby("raceid").apply(lambda xs: solve_ps(t, xs["winodds"]))
    ps = dt.groupby("raceid")["winodds"].apply(lambda xs: solve_ps(t, xs["winodds"]))
    ps.index = dt.index
    return np.sum(np.log(ps[dt.winner==1]))


def loglike2(t, grpd, winner_indx):
    ps = grpd.apply(lambda xs: solve_ps(t, xs))
    return np.sum(np.log(ps[winner_indx == 1]))

def main():
    winner_indx = df.winner == 1
    grpd = df.groupby("raceid")["winodds"]
    return loglike2(0.5, grpd, winner_indx)

%timeit print main()

foo = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                   'B' : ['one', 'one', 'two', 'three',
                          'two', 'two', 'one', 'three'],
                   'C' : np.random.randn(8),
                   'D' : np.random.randn(8)})

dg = df.groupby("A")

def foo(xs):
    return solve_ps(0.5, xs)

dg.apply(lambda x: np.sum(x))

