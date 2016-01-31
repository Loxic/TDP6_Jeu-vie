#TESTS sequentiels
N=30

echo "TESTS TAILLE ENTREE"

echo "TESTS SEQUENTIEL"
printf "" > test_seq.dat
for i in {1000..20000..1000}
do
    printf "$i " >> test_seq.dat
    salloc -N 1 ./life_seq $N $i >> test_seq.dat
    echo "test seq $i"
done

echo "TESTS OPENMP"
printf "" > test_omp.dat
for i in {1000..20000..1000}
do
    printf "$i " >> test_omp.dat
    salloc -N 1 ./life_omp $N $i >> test_omp.dat
    echo "test omp parts $i"
done

#echo "TESTS PTHREADS"
#printf "" > test_pthreads.dat
#for i in {1000..20000..1000}
#do
#    salloc -N 1 ./life_pthread $N $i >> test_seq_parts.dat
#    echo "test pthreads $i"
#done

echo "TESTS MPI SYNC"
printf "" > test_mpi_sync.dat
for i in {1000..20000..1000}
do
    printf "$i " >> test_mpi_sync.dat
    salloc -N 1 mpirun ./life_mpi_sync $N $i >> test_mpi_sync.dat
    echo "test mpi sync $i"
done

echo "TESTS MPI ASYNC"
printf "" > test_mpi_async.dat
for i in {1000..20000..1000}
do
    printf "$i " >> test_mpi_async.dat
    salloc -N 1 mpirun ./life_mpi_async $N $i >> test_mpi_async.dat
    echo "test mpi async $i"
done

echo "TESTS MPI PERSISTANT"
printf "" > test_mpi_pers.dat
for i in {1000..20000..1000}
do
    printf "$i " >> test_mpi_pers.dat
    salloc -N 1 mpirun ./life_mpi_pers $N $i >> test_mpi_pers.dat
    echo "test mpi_pers $i"
done


echo "TESTS SPEEDUP"
printf "" > test_speedup_seq.dat
salloc -N 1 ./life_seq $N 3000 >> test_speedup_seq.dat
echo "test speedup seq"


printf "" > test_speedup_mpi_sync.dat
for i in {1, 4, 9, 16, 25}
do
    printf "$i " >> test_speedup_mpi_sync.dat
    salloc -N 1 mpirun -n $i ./life_mpi_sync $N 3000 >> test_speedup_mpi_sync.dat
    echo "test speedup mpi_sync $i"
done


printf "" > test_speedup_mpi_async.dat
for i in {1, 4, 9, 16, 25}
do
    printf "$i " >> test_speed_mpi_async.dat
    salloc -N 1 mpirun -n $i ./life_mpi_sync $N 3000 >> test_speedup_mpi_async.dat
    echo "test speedup mpi_async $i"
done


printf "" > test_speedup_mpi_pers.dat
for i in {1, 4, 9, 16, 25}
do
    printf "$i " >> test_speedup_mpi_pers.dat
    salloc -N 1 mpirun -n $i ./life_mpi_sync $N 3000 >> test_speedup_mpi_pers.dat
    echo "test speedup mpi_pers $i"
done

printf "" > test_speedup_omp.dat
for i in {1, 4, 9, 16, 25}
do
    printf "$i " >> test_speedup_mpi_pers.dat
    salloc -N 1 ./life_omp $N 3000 $i >> test_speedup_mpi_pers.dat
    echo "test speedup omp $i"
done


