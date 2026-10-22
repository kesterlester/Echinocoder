# In general, data is a length-n mulitiset contianing "n" orderd lists, each of length-m.
# However, if m=1, the data can just be a length-n multiset of values.
# I.e. the following are OK:
#
# (m,n)=(1,4) : data=(1,2,3,4,)                  size=(4) or (1,4)
# (m,n)=(1,4) : data=((1),(2),(3),(4),)          size=(4,1)
# (m,n)=(2,4) : data=((1,2),(3,4),(5,6),(7,8) )  size=(4,2)


import numpy as np

def random_real_array_data(mn=(1,1) ):
    m,n = mn
    return np.random.randint(low=-9,high=10,size=(n,m))

# some random data that is a length-n list of reals
def random_real_linear_data(n=1):
    return np.random.randint(low=-9,high=10,size=n)

# some random data that is a length-n list of complex numbers
def random_complex_linear_data(n=1):
    return random_real_linear_data(n) + complex(0,1)*random_real_linear_data(n)




