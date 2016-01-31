#Plot SPEEDUP (3 courbes MPI, 1 openMP, 1 pthreads)
gnuplot <<EOF
set xlabel "Nombres de coeurs"
set ylabel "Temps d'execution (ms)"
set grid
set title "Speedup des différentes versions parallèles du jeu de la vie"
set terminal png
set output "speedup.png"
plot "test_speedup_mpi_sync.dat" title 'MPI SYNC' with lines, "test_speedup_mpi_async.dat" title 'MPI ASYNC' with lines, "test_speedup_mpi_pers.dat" title 'MPI PERSISTANT' with lines, "test_speedup_omp.dat" title 'OpenMP' with lines
EOF



#Plot comportement
gnuplot <<EOF
set xlabel "Nombre de cellules"
set ylabel "Temps d'execution (ms)"
set grid
set title "Courbe de temps de l'algorithme sequentiel en fonction de la taille de l'entrée"
set terminal png
set output "seq.png"
plot "test_seq.dat" title 'Sequentiel'
EOF


gnuplot <<EOF
set xlabel "Nombre de cellules"
set ylabel "Temps d'execution (ms)"
set grid
set title "Courbe de temps de l'algorithme OpenMP en fonction de la taille de l'entrée"
set terminal png
set output "omp.png"
plot "test_omp.dat" title 'OpenMP'
EOF

gnuplot <<EOF
set xlabel "Nombre de cellules"
set ylabel "Temps d'execution (ms)"
set grid
set title "Courbe de temps de l'algorithme MPI synchronisé en fonction de la taille de l'entrée"
set terminal png
set output "mpi_sync.png"
plot "test_mpi_sync.dat" title 'MPI SYNC'
EOF

gnuplot <<EOF
set xlabel "Nombre de cellules"
set ylabel "Temps d'execution (ms)"
set grid
set title "Courbe de temps de l'algorithme MPI asynchrone en fonction de la taille de l'entrée"
set terminal png
set output "mpi_async.png"
plot "test_mpi_async.dat" title 'MPI ASYNC'
EOF

gnuplot <<EOF
set xlabel "Nombre de cellules"
set ylabel "Temps d'execution (ms)"
set grid
set title "Courbe de temps de l'algorithme MPI persistant en fonction de la taille de l'entrée"
set terminal png
set output "mpi_pers.png"
plot "test_mpi_pers.dat" title 'MPI PERS'
EOF


