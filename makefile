f90 = gfortran
objects =2Deuler.o cellnode.o areacalc.o edgenormals.o edgearrange.o Roe_interior.o Roe_wallflux.o timestep.o
xfvm: $(objects)
	$(f90) $(objects)
	mv a.out xfvm
2Deuler.o : 2Deuler.f90
	$(f90) -c 2Deuler.f90
cellnode.o:	cellnode.f90
	$(f90) -c cellnode.f90
areacalc.o:	areacalc.f90
	$(f90) -c areacalc.f90
edgenormals.o:	edgenormals.f90
	$(f90) -c edgenormals.f90
edgearrange.o: edgearrange.f90
	$(f90) -c edgearrange.f90
Roe_interior.o: Roe_interior.f90
	$(f90) -c Roe_interior.f90
Roe_wallflux.o: Roe_wallflux.f90
	$(f90) -c Roe_wallflux.f90
timestep.o: timestep.f90
	$(f90) -c timestep.f90
clean:
	rm $(objects)

