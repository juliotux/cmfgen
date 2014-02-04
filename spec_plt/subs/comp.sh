#
set FFLAG=(-w -Fixed -X 9 -Am)
#
f90 -c $FFLAG convolve.f
f90 -c $FFLAG h2abs.f
f90 -c $FFLAG hiabs.f
f90 -c $FFLAG num_rec.f
f90 -c $FFLAG uvabs.f
#
ar -ruv ../../lib/libspec.a *.o
