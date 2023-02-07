TEST: main.cpp intersections.h intersections.cpp fastbbf.h fastbbf.cpp util.h util.cpp upcStatistical.h upcStatistical.cpp
	g++ -O3 -fopenmp main.cpp intersections.h intersections.cpp fastbbf.h fastbbf.cpp util.h util.cpp upcStatistical.h upcStatistical.cpp -o TEST
run:	TEST
	./TEST