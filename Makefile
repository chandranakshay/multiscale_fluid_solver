CC = mpicxx 
CFLAGS =-c -O3 -std=c++11 -mavx -DAVX -I . -I . -I ./RD3Q41 -I ./moldiceRNG/include -I ./src -I ./RD3Q41/parallelisation -I ./RD3Q41/parallelisation/noiseCommunication  -I ./RD3Q41/parallelisation/partialCommunication -I ./DSMC_Solver
all : haha 

haha : main.o procInfo3D.o D3RQ41.o pack41.o partial_pack41.o marker.o equilibrium41.o initialConditions.o collide41.o ioStreamGrid.o communication41.o partial_communication41.o boundaryCondition.o inletOutlet.o advection41.o bounceBack41_withForce.o diffuseBounceBack41_withForce.o diffuse41.o pvts.o restart.o momentsForDSMC.o perlinNoise.o noiseFunctions.o
	$(CC) main.o procInfo3D.o D3RQ41.o pack41.o partial_pack41.o marker.o equilibrium41.o initialConditions.o collide41.o ioStreamGrid.o communication41.o partial_communication41.o boundaryCondition.o inletOutlet.o advection41.o bounceBack41_withForce.o diffuseBounceBack41_withForce.o diffuse41.o pvts.o restart.o momentsForDSMC.o perlinNoise.o noiseFunctions.o -o baseCoupled.out 

main.o : ./DSMC_Solver/main.C
	$(CC) $(CFLAGS) ./DSMC_Solver/main.C 

procInfo3D.o : ./src/procInfo3D.C
	$(CC) $(CFLAGS) ./src/procInfo3D.C

D3RQ41.o : ./src/D3RQ41.C
	$(CC) $(CFLAGS) ./src/D3RQ41.C

pack41.o : ./RD3Q41/parallelisation/pack41.C
	$(CC) $(CFLAGS) ./RD3Q41/parallelisation/pack41.C

partial_pack41.o : ./RD3Q41/parallelisation/partialCommunication/partial_pack41.C
	$(CC) $(CFLAGS) ./RD3Q41/parallelisation/partialCommunication/partial_pack41.C

marker.o : ./RD3Q41/marker.C
	$(CC) $(CFLAGS) ./RD3Q41/marker.C 

equilibrium41.o : ./RD3Q41/equilibrium41.C
	$(CC) $(CFLAGS) ./RD3Q41/equilibrium41.C

initialConditions.o : ./RD3Q41/initialConditions.C
	$(CC) $(CFLAGS) ./RD3Q41/initialConditions.C

collide41.o : ./RD3Q41/collide41.C
	$(CC) $(CFLAGS) ./RD3Q41/collide41.C

communication41.o : ./RD3Q41/parallelisation/communication41.C
	$(CC) $(CFLAGS) ./RD3Q41/parallelisation/communication41.C

partial_communication41.o : ./RD3Q41/parallelisation/partialCommunication/partial_communication41.C
	$(CC) $(CFLAGS) ./RD3Q41/parallelisation/partialCommunication/partial_communication41.C

ioStreamGrid.o : ./RD3Q41/parallelisation/ioStreamGrid.C
	$(CC) $(CFLAGS) ./RD3Q41/parallelisation/ioStreamGrid.C

boundaryCondition.o : ./RD3Q41/boundaryCondition.C
	$(CC) $(CFLAGS) ./RD3Q41/boundaryCondition.C

inletOutlet.o : ./RD3Q41/inletOutlet.C
	$(CC) $(CFLAGS) ./RD3Q41/inletOutlet.C

advection41.o : ./RD3Q41/advection41.C
	$(CC) $(CFLAGS) ./RD3Q41/advection41.C

bounceBack41_withForce.o : ./RD3Q41/bounceBack41_withForce.C
	$(CC) $(CFLAGS) ./RD3Q41/bounceBack41_withForce.C

diffuseBounceBack41_withForce.o : ./RD3Q41/diffuseBounceBack41_withForce.C
	$(CC) $(CFLAGS) ./RD3Q41/diffuseBounceBack41_withForce.C

	
diffuse41.o : ./RD3Q41/diffuse41.C
	$(CC) $(CFLAGS) ./RD3Q41/diffuse41.C

pvts.o : ./RD3Q41/pvts.C
	$(CC) $(CFLAGS) ./RD3Q41/pvts.C

restart.o : ./RD3Q41/restart.C
	$(CC) $(CFLAGS) ./RD3Q41/restart.C

momentsForDSMC.o : ./RD3Q41/momentsForDSMC.C
	$(CC) $(CFLAGS) ./RD3Q41/momentsForDSMC.C

# bounceBack41.o : ./RD3Q41/bounceBack41.C
# 	$(CC) $(CFLAGS) ./RD3Q41/bounceBack41.C
# 

perlinNoise.o : ./RD3Q41/perlinNoise.C
	$(CC) $(CFLAGS) ./RD3Q41/perlinNoise.C

noiseFunctions.o : ./RD3Q41/noiseFunctions.C
	$(CC) $(CFLAGS) ./RD3Q41/noiseFunctions.C
 	
	
clean:
	rm -rf mass*.txt enstrophy.txt output.dat *.out *.o results/* *.err *.out *.txt moments/* ./momentsLocal/* ./momentsFlux/* ./symbols ./c_test* ./core*
