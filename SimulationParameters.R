########################################################
########################################################

#Data generating process setting
method =  'Jose'
#method = 'McG'

#number of covariates
p = 50

#Number of variables that conforms the denser network

T0 = 30

#number of samples
n = 1000

#Number of hub nodes
r = 3

#probability of connection for a hub node
ph = 0.9

#probability of connection between non-hub nodes
pnh = 0.05

#probability of connection for nodes outside the T0xT0 matrix
pneff = 0.001

shuffle = FALSE
# a, b range for the values of the precision matrix
hmin = 4
hmax = 6

nhmin = 4
nhmax = 6

neffmin = 4
neffmax = 6

# Delta, minimum eigenvalue separation from zero
diagonal_shift = 5

seed = 12

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#######################################################################

args = list(method = method, p = p, T0 = T0, n = n, r = r, ph = ph, pnh = pnh, pneff = pneff, shuffle = shuffle,
            hmin = hmin, hmax = hmax, nhmin = nhmin, nhmax = nhmax, neffmin = neffmin, neffmax = neffmax, diagonal_shift = diagonal_shift, seed = seed)
