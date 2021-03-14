# Meeting : 03-07-21
# Create a matrix with all combinations of Xbin, Ybin, Tbin
# Discussed taking Croneker product of the 3 matricies

#Practice from meeting
Ndays = 5 # Fill this in based on computation in MRF code
Nbins = 3

N = Ndays*Nbins*Nbins
N

T = rep(0,N)
X=T
Y=T
data.frame(T,X,Y)

for (t in 1:Ndays){
  T[((t-1)*Nbins*Nbins + 1) : (t*Nbins*Nbins)]  = t
  for (x in 1:Nbins) {
    X[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = x #Give all Values for X
    Y[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = 1:Nbins #Give all Values for Y
  }
}

Z = rep(0,N)
Z[(t-1)*Nbins^2 + (x-1)*Nbins + y] = 1
data.frame(T,X,Y,Z)

# Test : Line number to check for event
t=2
x=2
y=2
(t-1)*Nbins^2 + (x-1)*Nbins + y #results on line 14
